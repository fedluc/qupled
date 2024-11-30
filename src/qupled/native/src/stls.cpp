#include "stls.hpp"
#include "bin_util.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "vector_util.hpp"
#include <fmt/core.h>

using namespace std;
using namespace vecUtil;
using namespace binUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;
using ItgType = Integrator1D::Type;

// -----------------------------------------------------------------
// STLS class
// -----------------------------------------------------------------

Stls::Stls(const StlsInput &in_, const bool verbose_, const bool writeFiles_)
    : Rpa(in_, verbose_),
      in(in_),
      writeFiles(writeFiles_ && isRoot()) {
  // Check if iet scheme should be solved
  useIet = in.getTheory() == "STLS-HNC" || in.getTheory() == "STLS-IOI"
           || in.getTheory() == "STLS-LCT";
  // Set name of recovery files
  recoveryFileName = fmt::format("recovery_rs{:.3f}_theta{:.3f}_{}.bin",
                                 in.getCoupling(),
                                 in.getDegeneracy(),
                                 in.getTheory());
  // Allocate arrays
  const size_t nx = wvg.size();
  slfcNew.resize(nx);
  if (useIet) { bf.resize(nx); }
}

int Stls::compute() {
  try {
    init();
    println("Structural properties calculation ...");
    doIterations();
    println("Done");
    return 0;
  } catch (const runtime_error &err) {
    cerr << err.what() << endl;
    return 1;
  }
}

// Initialize basic properties
void Stls::init() {
  Rpa::init();
  if (useIet) {
    print("Computing bridge function adder: ");
    computeBf();
    println("Done");
  }
}

// Compute static local field correction
void Stls::computeSlfc() {
  assert(ssf.size() == wvg.size());
  assert(slfc.size() == wvg.size());
  computeSlfcStls();
  if (useIet) computeSlfcIet();
}

void Stls::computeSlfcStls() {
  const int nx = wvg.size();
  const Interpolator1D itp(wvg, ssf);
  for (int i = 0; i < nx; ++i) {
    Slfc slfcTmp(wvg[i], wvg.front(), wvg.back(), itp, itg);
    slfcNew[i] = slfcTmp.get();
  }
}

void Stls::computeSlfcIet() {
  Integrator2D itg2(in.getIntError());
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  const Interpolator1D ssfItp(wvg, ssf);
  const Interpolator1D slfcItp(wvg, slfc);
  const Interpolator1D bfItp(wvg, bf);
  for (size_t i = 0; i < wvg.size(); ++i) {
    SlfcIet slfcTmp(
        wvg[i], wvg.front(), wvg.back(), ssfItp, slfcItp, bfItp, itgGrid, itg2);
    slfcNew[i] += slfcTmp.get();
  }
}

// Compute bridge function
void Stls::computeBf() {
  const size_t nx = wvg.size();
  Integrator1D itgF(ItgType::FOURIER, 1e-10);
  assert(bf.size() == nx);
  for (size_t i = 0; i < nx; ++i) {
    BridgeFunction bfTmp(in.getTheory(),
                         in.getIETMapping(),
                         in.getCoupling(),
                         in.getDegeneracy(),
                         wvg[i],
                         itgF);
    bf[i] = bfTmp.get();
  }
}

// stls iterations
void Stls::doIterations() {
  const int maxIter = in.getNIter();
  const int outIter = in.getOutIter();
  const double minErr = in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  initialGuess();
  while (counter < maxIter + 1 && err > minErr) {
    // Start timing
    double tic = timer();
    // Update static structure factor
    computeSsf();
    // Update static local field correction
    computeSlfc();
    // Update diagnostic
    counter++;
    err = computeError();
    // Update solution
    updateSolution();
    // Write output
    if (counter % outIter == 0 && writeFiles) { writeRecovery(); }
    // End timing
    double toc = timer();
    // Print diagnostic
    println(fmt::format("--- iteration {:d} ---", counter));
    println(fmt::format("Elapsed time: {:.3f} seconds", toc - tic));
    println(fmt::format("Residual error: {:.5e}", err));
    fflush(stdout);
  }
}

// Initial guess for stls iterations
void Stls::initialGuess() {
  // From recovery file
  if (initialGuessFromRecovery()) { return; }
  // From guess in input
  if (initialGuessFromInput()) { return; }
  // Default
  fill(slfc.begin(), slfc.end(), 0.0);
}

bool Stls::initialGuessFromRecovery() {
  vector<double> wvgFile;
  vector<double> slfcFile;
  readRecovery(wvgFile, slfcFile);
  const Interpolator1D slfci(wvgFile, slfcFile);
  if (!slfci.isValid()) { return false; }
  const double xmaxi = wvgFile.back();
  for (size_t i = 0; i < wvg.size(); ++i) {
    const double x = wvg[i];
    if (x <= xmaxi) {
      slfc[i] = slfci.eval(x);
    } else {
      slfc[i] = 1.0;
    }
  }
  return true;
}

bool Stls::initialGuessFromInput() {
  const Interpolator1D slfci(in.getGuess().wvg, in.getGuess().slfc);
  if (!slfci.isValid()) { return false; }
  const double xmaxi = in.getGuess().wvg.back();
  for (size_t i = 0; i < wvg.size(); ++i) {
    const double x = wvg[i];
    if (x <= xmaxi) {
      slfc[i] = slfci.eval(x);
    } else {
      slfc[i] = 1.0;
    }
  }
  return true;
}

// Compute residual error for the stls iterations
double Stls::computeError() const { return rms(slfc, slfcNew, false); }

// Update solution during stls iterations
void Stls::updateSolution() {
  const double aMix = in.getMixingParameter();
  slfc = linearCombination(slfcNew, aMix, slfc, 1 - aMix);
}

// Recovery files
void Stls::writeRecovery() {
  ofstream file;
  file.open(recoveryFileName, ios::binary);
  if (!file.is_open()) {
    throwError("Recovery file " + recoveryFileName + " could not be created.");
  }
  int nx = wvg.size();
  writeDataToBinary<int>(file, nx);
  writeDataToBinary<decltype(wvg)>(file, wvg);
  writeDataToBinary<decltype(slfc)>(file, slfc);
  file.close();
  if (!file) {
    throwError("Error in writing the recovery file " + recoveryFileName);
  }
}

void Stls::readRecovery(vector<double> &wvgFile,
                        vector<double> &slfcFile) const {
  const string fileName = in.getRecoveryFileName();
  if (fileName.empty()) { return; }
  ifstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throwError("Output file " + fileName + " could not be opened.");
  }
  int nx;
  readDataFromBinary<int>(file, nx);
  wvgFile.resize(nx);
  slfcFile.resize(nx);
  readDataFromBinary<decltype(wvgFile)>(file, wvgFile);
  readDataFromBinary<decltype(slfcFile)>(file, slfcFile);
  file.close();
  if (!file) { throwError("Error in reading from file " + fileName); }
}

// -----------------------------------------------------------------
// SlfcBase class
// -----------------------------------------------------------------

// Compute static structure factor from interpolator
double SlfcBase::ssf(const double &y) const { return ssfi.eval(y); }

// -----------------------------------------------------------------
// Slfc class
// -----------------------------------------------------------------

// Get result of integration
double Slfc::get() const {
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg.compute(func, ItgParam(yMin, yMax));
  return itg.getSolution();
}

// Integrand
double Slfc::integrand(const double &y) const {
  double y2 = y * y;
  double x2 = x * x;
  if (x == 0.0 || y == 0.0) { return 0.0; }
  if (x == y) { return -(3.0 / 4.0) * y2 * (ssf(y) - 1.0); };
  if (x > y) {
    return -(3.0 / 4.0) * y2 * (ssf(y) - 1.0)
           * (1 + (x2 - y2) / (2 * x * y) * log((x + y) / (x - y)));
  }
  return -(3.0 / 4.0) * y2 * (ssf(y) - 1.0)
         * (1 + (x2 - y2) / (2 * x * y) * log((x + y) / (y - x)));
}

// -----------------------------------------------------------------
// SlfcIet class
// -----------------------------------------------------------------

// Compute static local field correction from interpolator
double SlfcIet::slfc(const double &y) const { return slfci.eval(y); }

// Compute bridge function from interpolator
double SlfcIet::bf(const double &y) const { return bfi.eval(y); }

// Get at finite temperature
double SlfcIet::get() const {
  if (x == 0.0) return 0.0;
  auto wMin = [&](const double &y) -> double {
    return (y > x) ? y - x : x - y;
  };
  auto wMax = [&](const double &y) -> double { return min(yMax, x + y); };
  auto func1 = [&](const double &y) -> double { return integrand1(y); };
  auto func2 = [&](const double &w) -> double { return integrand2(w); };
  itg.compute(func1, func2, Itg2DParam(yMin, yMax, wMin, wMax), itgGrid);
  return 3.0 / (8.0 * x) * itg.getSolution() + bf(x);
}

// Level 1 integrand
double SlfcIet::integrand1(const double &y) const {
  if (y == 0.0) return 0.0;
  return (-bf(y) - (ssf(y) - 1.0) * (slfc(y) - 1.0)) / y;
}

// Level 2 integrand
double SlfcIet::integrand2(const double &w) const {
  const double y = itg.getX();
  const double y2 = y * y;
  const double w2 = w * w;
  const double x2 = x * x;
  return (w2 - y2 - x2) * w * (ssf(w) - 1.0);
}

// -----------------------------------------------------------------
// BridgeFunction class
// -----------------------------------------------------------------

double BridgeFunction::get() const {
  if (theory == "STLS-HNC" || theory == "QSTLS-HNC") { return hnc(); }
  if (theory == "STLS-IOI" || theory == "QSTLS-IOI") { return ioi(); }
  if (theory == "STLS-LCT" || theory == "QSTLS-LCT") { return lct(); }
  throwError("Unknown theory to compute the bridge function term");
  return numUtil::Inf;
}

double BridgeFunction::couplingParameter() const {
  const double fact = 2 * lambda * lambda * rs;
  if (mapping == "sqrt") { return fact / sqrt(1 + Theta * Theta); }
  if (mapping == "linear") { return fact / (1 + Theta); }
  if (Theta != 0.0) { return fact / Theta; }
  throwError("The standard iet mapping cannot be used in the "
             "ground state");
  return numUtil::Inf;
}

double BridgeFunction::hnc() const { return 0.0; }

double BridgeFunction::ioi() const {
  const double l2 = lambda * lambda;
  const double l3 = l2 * lambda;
  const double l4 = l3 * lambda;
  const double l5 = l4 * lambda;
  const double l6 = l5 * lambda;
  const double l7 = l6 * lambda;
  const double l8 = l7 * lambda;
  const double Gamma = couplingParameter();
  const double lnG = log(Gamma);
  const double lnG2 = lnG * lnG;
  const double b0 = 0.258 - 0.0612 * lnG + 0.0123 * lnG2 - 1.0 / Gamma;
  const double b1 = 0.0269 + 0.0318 * lnG + 0.00814 * lnG2;
  if (b0 / b1 <= 0.0 || Gamma < 5.25 || Gamma > 171.8) {
    const string msg = fmt::format("The IET schemes cannot be applied "
                                   "to this state point because Gamma = {:.8f} "
                                   "falls outside the range of validty of the "
                                   "bridge function parameterization\n",
                                   Gamma);
    throwError(msg);
  }
  const double c1 = 0.498 - 0.280 * lnG + 0.0294 * lnG2;
  const double c2 = -0.412 + 0.219 * lnG - 0.0251 * lnG2;
  const double c3 = 0.0988 - 0.0534 * lnG + 0.00682 * lnG2;
  const double b02 = b0 * b0;
  const double b03 = b02 * b0;
  const double b04 = b03 * b0;
  const double b05 = b04 * b0;
  const double b06 = b05 * b0;
  const double b07 = b06 * b0;
  const double b08 = b07 * b0;
  const double b12 = b1 * b1;
  const double b13 = b12 * b1;
  const double b14 = b13 * b1;
  const double b15 = b14 * b1;
  const double b16 = b15 * b1;
  const double b17 = b16 * b1;
  const double b18 = b17 * b1;
  const double b02_b12 = b02 / b12;
  const double b03_b13 = b03 / b13;
  const double b04_b14 = b04 / b14;
  const double b05_b15 = b05 / b15;
  const double b06_b16 = b06 / b16;
  const double b07_b17 = b07 / b17;
  const double b08_b18 = b08 / b18;
  const double fact = sqrt(M_PI) / (4.0 * l2) * pow(b0 / b1, 1.5);
  const double q2 = x * x;
  const double q3 = q2 * x;
  const double q4 = q3 * x;
  const double q5 = q4 * x;
  const double q6 = q5 * x;
  const double q7 = q6 * x;
  const double q8 = q7 * x;
  const double bf1 =
      -b0
      + c1 / 16.0
            * (60.0 * b02_b12 - 20.0 * b03_b13 * q2 / l2 + b04_b14 * q4 / l4);
  const double bf2 = c2 / 64.0
                     * (840.0 * b03_b13 - 420.0 * b04_b14 * q2 / l2
                        + 42.0 * b05_b15 * q4 / l4 - b06_b16 * q6 / l6);
  ;
  const double bf3 = c3 / 256.0
                     * (15120.0 * b04_b14 - 10080.0 * b05_b15 * q2 / l2
                        + 1512.0 * b06_b16 * q4 / l4 - 72.0 * b07_b17 * q6 / l6
                        + b08_b18 * q8 / l8);
  return fact * q2 * (bf1 + bf2 + bf3) * exp(-b0 * q2 / (4.0 * b1 * l2));
}

double BridgeFunction::lct() const {
  const double Gamma = couplingParameter();
  auto func = [&](const double &r) -> double { return lctIntegrand(r, Gamma); };
  itg.compute(func, ItgParam(x / lambda));
  return itg.getSolution() * (x / lambda) / Gamma;
  return 0.0;
}

double BridgeFunction::lctIntegrand(const double &r,
                                    const double &Gamma) const {
  if (Gamma < 5.0) {
    const string msg = fmt::format("The IET schemes cannot be applied "
                                   "to this state point because Gamma = {:.8f} "
                                   "falls outside the range of validty of the "
                                   "bridge function parameterization\n",
                                   Gamma);
    throwError(msg);
  }
  const double Gamma1_6 = pow(Gamma, 1. / 6.);
  const double lnG = log(Gamma);
  const double lnG2 = lnG * lnG;
  const double lnG3 = lnG2 * lnG;
  const double lnG4 = lnG3 * lnG;
  const double a0 =
      Gamma * (0.076912 - 0.10465 * lnG + 0.0056629 * lnG2 + 0.00025656 * lnG3);
  const double a2 =
      Gamma * (0.068045 - 0.036952 * lnG + 0.048818 * lnG2 - 0.0048985 * lnG3);
  const double a3 =
      Gamma * (-0.30231 + 0.30457 * lnG - 0.11424 * lnG2 + 0.0095993 * lnG3);
  const double a4 =
      Gamma * (0.25111 - 0.26800 * lnG + 0.082268 * lnG2 - 0.0064960 * lnG3);
  const double a5 =
      Gamma * (-0.061894 + 0.066811 * lnG - 0.019140 * lnG2 + 0.0014743 * lnG3);
  const double c0 = Gamma
                    * (0.25264 - 0.31615 * lnG + 0.13135 * lnG2
                       - 0.023044 * lnG3 + 0.0014666 * lnG4);
  const double c1 = Gamma1_6
                    * (-12.665 + 20.802 * lnG - 9.6296 * lnG2 + 1.7889 * lnG3
                       - 0.11810 * lnG4);
  const double c2 = Gamma1_6
                    * (15.285 - 14.076 * lnG + 5.7558 * lnG2 - 1.0188 * lnG3
                       + 0.06551 * lnG4);
  const double c3 = Gamma1_6
                    * (35.330 - 40.727 * lnG + 16.690 * lnG2 - 2.8905 * lnG3
                       + 0.18243 * lnG4);
  const double r2 = r * r;
  const double r3 = r2 * r;
  const double r4 = r3 * r;
  const double r5 = r4 * r;
  const double rshift = r - 1.44;
  const double bsr = a0 + a2 * r2 + a3 * r3 + a4 * r4 + a5 * r5;
  const double blr = c0 * exp(-c1 * rshift) * exp(-0.3 * r2)
                     * (cos(c2 * rshift) + c3 * exp(-3.5 * rshift));
  const double sf = 0.5 * (1.0 + erf(5.0 * (r - 1.50)));
  return r * ((1 - sf) * bsr + sf * blr);
}
