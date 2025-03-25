#include "stls2d.hpp"
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

Stls2D::Stls2D(const StlsInput &in_, const bool verbose_, const bool writeFiles_)
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
