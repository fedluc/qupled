#include "qstls.hpp"
#include "bin_util.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "vector_util.hpp"
#include <filesystem>
#include <fmt/core.h>
#include <numeric>

using namespace std;
using namespace vecUtil;
using namespace binUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;

// -----------------------------------------------------------------
// QSTLS class
// -----------------------------------------------------------------

Qstls::Qstls(const QstlsInput &in_, const bool verbose_, const bool writeFiles_)
    : Stls(in_, verbose_, writeFiles_),
      in(in_) {
  // Throw error message for ground state calculations
  if (in.getDegeneracy() == 0.0 && useIet) {
    throwError("Ground state calculations are not available "
               "for the quantum IET schemes");
  }
  // Check if iet scheme should be solved
  useIet = in.getTheory() == "QSTLS-HNC" || in.getTheory() == "QSTLS-IOI"
           || in.getTheory() == "QSTLS-LCT";
  // Allocate arrays
  const size_t nx = wvg.size();
  const size_t nl = in.getNMatsubara();
  adr.resize(nx, nl);
  ssfNew.resize(nx);
  ssfOld.resize(nx);
  adrFixed.resize(nx, nl, nx);
  if (useIet) {
    bf.resize(nx);
    adrOld.resize(nx, nl);
  }
  // Deallocate arrays that are inherited but not used
  vector<double>().swap(slfcNew);
}

int Qstls::compute() {
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

void Qstls::init() {
  Stls::init();
  print("Computing fixed component of the auxiliary density response: ");
  computeAdrFixed();
  println("Done");
  if (useIet) {
    print("Computing fixed component of the iet auxiliary density response: ");
    computeAdrFixedIet();
    println("Done");
  }
}

// qstls iterations
void Qstls::doIterations() {
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
    // Update auxiliary density response
    computeAdr();
    // Update static structure factor
    computeSsf();
    // Update diagnostic
    counter++;
    err = computeError();
    // Update solution
    updateSolution();
    // Write output
    if (counter % outIter == 0 && writeFiles) { writeRecovery(); };
    // End timing
    double toc = timer();
    // Print diagnostic
    println(fmt::format("--- iteration {:d} ---", counter));
    println(fmt::format("Elapsed time: {:.3f} seconds", toc - tic));
    println(fmt::format("Residual error: {:.5e}", err));
    fflush(stdout);
  }
  // Set static structure factor for output
  ssf = ssfOld;
}

// Initial guess for qstls iterations
void Qstls::initialGuess() {
  assert(!ssfNew.empty());
  assert(!ssfOld.empty());
  assert(!adr.empty());
  assert(!useIet || !adrOld.empty());
  // From recovery file
  if (initialGuessFromRecovery()) { return; }
  // From guess in input
  if (initialGuessFromInput()) { return; }
  // Default
  Rpa rpa(in, false);
  int status = rpa.compute();
  if (status != 0) {
    throwError("Failed to compute the default initial guess");
  }
  ssfOld = rpa.getSsf();
  if (useIet) { adrOld.fill(0.0); }
}

bool Qstls::initialGuessFromRecovery() {
  vector<double> wvg_;
  vector<double> ssf_;
  Vector2D adr_;
  Vector3D adrFixed_;
  double Theta;
  int nl_;
  readRecovery(wvg_, ssf_, adr_, adrFixed_, Theta, nl_);
  const bool ssfIsSet = initialGuessSsf(wvg_, ssf_);
  const bool adrIsSet = (useIet) ? initialGuessAdr(wvg_, adr_) : true;
  const bool adrFixedIsSet = initialGuessAdrFixed(wvg_, Theta, nl_, adrFixed_);
  return ssfIsSet && adrIsSet && adrFixedIsSet;
}

bool Qstls::initialGuessFromInput() {
  const auto &guess = in.getGuess();
  const bool ssfIsSet = initialGuessSsf(guess.wvg, guess.ssf);
  const bool adrIsSet = (useIet) ? initialGuessAdr(guess.wvg, guess.adr) : true;
  return ssfIsSet && adrIsSet;
}

bool Qstls::initialGuessSsf(const vector<double> &wvg_,
                            const vector<double> &ssf_) {
  const int nx = wvg.size();
  const Interpolator1D ssfi(wvg_, ssf_);
  if (!ssfi.isValid()) { return false; }
  const double xMax = wvg_.back();
  for (int i = 0; i < nx; ++i) {
    const double &x = wvg[i];
    ssfOld[i] = (x <= xMax) ? ssfi.eval(x) : 1.0;
  }
  return true;
}

bool Qstls::initialGuessAdr(const vector<double> &wvg_, const Vector2D &adr_) {
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const int nx_ = adr_.size(0);
  const int nl_ = adr_.size(1);
  const double xMax = (wvg_.empty()) ? 0.0 : wvg_.back();
  vector<Interpolator1D> itp(nl_);
  for (int l = 0; l < nl_; ++l) {
    vector<double> tmp(nx_);
    for (int i = 0; i < nx_; ++i) {
      tmp[i] = adr_(i, l);
    }
    itp[l].reset(wvg_[0], tmp[0], nx_);
    if (!itp[l].isValid()) { return false; }
  }
  for (int i = 0; i < nx; ++i) {
    const double &x = wvg[i];
    if (x > xMax) {
      adrOld.fill(i, 0.0);
      continue;
    }
    for (int l = 0; l < nl; ++l) {
      adrOld(i, l) = (l < nl_) ? itp[l].eval(x) : 0.0;
    }
  }
  return true;
}

bool Qstls::initialGuessAdrFixed(const vector<double> &wvg_,
                                 const double &Theta,
                                 const int &nl_,
                                 const Vector3D &adrFixed_) {
  if (!checkAdrFixed(wvg_, Theta, nl_)) { return false; }
  adrFixed = adrFixed_;
  return true;
}

// Compute auxiliary density response
void Qstls::computeAdr() {
  if (in.getDegeneracy() == 0.0) { return; }
  const int nx = wvg.size();
  const Interpolator1D ssfi(wvg, ssfOld);
  for (int i = 0; i < nx; ++i) {
    Adr adrTmp(in.getDegeneracy(), wvg.front(), wvg.back(), wvg[i], ssfi, itg);
    adrTmp.get(wvg, adrFixed, adr);
  }
  if (useIet) computeAdrIet();
  for (int i = 0; i < nx; ++i) {
    slfc[i] = adr(i, 0) / idr(i, 0);
  };
}

// Compute static structure factor
void Qstls::computeSsf() {
  if (in.getDegeneracy() == 0.0) {
    computeSsfGround();
  } else {
    computeSsfFinite();
  }
}

// Compute static structure factor at finite temperature
void Qstls::computeSsfFinite() {
  if (in.getCoupling() > 0.0) {
    assert(ssfNew.size() > 0);
    assert(adr.size() > 0);
    assert(idr.size() > 0);
    if (useIet) assert(bf.size() > 0);
  }
  const double Theta = in.getDegeneracy();
  const double rs = in.getCoupling();
  const int nx = wvg.size();
  const int nl = idr.size(1);
  for (int i = 0; i < nx; ++i) {
    const double bfi = (useIet) ? bf[i] : 0;
    Qssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], nl, &idr(i), &adr(i), bfi);
    ssfNew[i] = ssfTmp.get();
  }
}

// Compute static structure factor at zero temperature
void Qstls::computeSsfGround() {
  const Interpolator1D ssfi(wvg, ssfOld);
  const double rs = in.getCoupling();
  const double OmegaMax = in.getFrequencyCutoff();
  const size_t nx = wvg.size();
  const double xMax = wvg.back();
  auto loopFunc = [&](int i) -> void {
    Integrator1D itgTmp(itg);
    QssfGround ssfTmp(wvg[i], rs, ssfHF[i], xMax, OmegaMax, ssfi, itgTmp);
    ssfNew[i] = ssfTmp.get();
  };
  const auto &loopData = parallelFor(loopFunc, nx, in.getNThreads());
  gatherLoopData(ssfNew.data(), loopData, nx);
}

// Compute residual error for the qstls iterations
double Qstls::computeError() const { return rms(ssfNew, ssfOld, false); }

// Update solution during qstls iterations
void Qstls::updateSolution() {
  const double aMix = in.getMixingParameter();
  ssfOld = linearCombination(ssfNew, aMix, ssfOld, 1 - aMix);
  if (useIet) {
    adrOld.mult(1 - aMix);
    adrOld.linearCombination(adr, aMix);
  }
}

void Qstls::computeAdrFixed() {
  if (in.getDegeneracy() == 0.0) { return; }
  // Check if it adrFixed can be loaded from input
  if (!in.getFixed().empty()) {
    readAdrFixedFile(adrFixed, in.getFixed(), false);
    return;
  }
  // Compute from scratch
  fflush(stdout);
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const int nxnl = nx * nl;
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  // Parallel for loop (Hybrid MPI and OpenMP)
  auto loopFunc = [&](int i) -> void {
    Integrator2D itg2(in.getIntError());
    AdrFixed adrTmp(
        in.getDegeneracy(), wvg.front(), wvg.back(), wvg[i], mu, itgGrid, itg2);
    adrTmp.get(wvg, adrFixed);
  };
  const auto &loopData = parallelFor(loopFunc, nx, in.getNThreads());
  gatherLoopData(adrFixed.data(), loopData, nxnl);
  // Write result to output file
  if (isRoot()) {
    try {
      writeAdrFixedFile(adrFixed, adrFixedFileName);
    } catch (...) {
      throwError("Error in the output file for the fixed component"
                 " of the auxiliary density response.");
    }
  }
}

void Qstls::writeAdrFixedFile(const Vector3D &res,
                              const string &fileName) const {
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const double Theta = in.getDegeneracy();
  ofstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throwError("Output file " + fileName + " could not be created.");
  }
  writeDataToBinary<int>(file, nx);
  writeDataToBinary<int>(file, nl);
  writeDataToBinary<double>(file, Theta);
  writeDataToBinary<vector<double>>(file, wvg);
  writeDataToBinary<Vector3D>(file, res);
  file.close();
  if (!file) { throwError("Error in writing to file " + fileName); }
}

void Qstls::readAdrFixedFile(Vector3D &res,
                             const string &fileName,
                             const bool iet) const {
  if (fileName.empty()) { return; }
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  ifstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throwError("Input file " + fileName + " could not be opened.");
  }
  int nx_;
  int nl_;
  vector<double> wvg_;
  double Theta_;
  readDataFromBinary<int>(file, nx_);
  readDataFromBinary<int>(file, nl_);
  readDataFromBinary<double>(file, Theta_);
  wvg_.resize(nx);
  readDataFromBinary<vector<double>>(file, wvg_);
  if (iet) {
    res.resize(nl, nx, nx);
  } else {
    res.resize(nx, nl, nx);
  }
  readDataFromBinary<Vector3D>(file, res);
  file.close();
  if (!file) { throwError("Error in reading from file " + fileName); }
  if (!checkAdrFixed(wvg_, Theta_, nl_)) {
    throwError("Fixed component of the auxiliary density response"
               " loaded from file is incompatible with input");
  }
}

bool Qstls::checkAdrFixed(const vector<double> &wvg_,
                          const double Theta_,
                          const int nl_) const {
  constexpr double tol = 1e-15;
  bool consistentGrid = wvg_.size() == wvg.size();
  if (consistentGrid) {
    const vector<double> wvgDiff = diff(wvg_, wvg);
    const double &wvgMaxDiff =
        abs(*max_element(wvgDiff.begin(), wvgDiff.end()));
    consistentGrid = consistentGrid && wvgMaxDiff <= tol;
  }
  const bool consistentMatsubara = nl_ == in.getNMatsubara();
  const bool consistentTheta = abs(Theta_ - in.getDegeneracy()) <= tol;
  return consistentMatsubara && consistentTheta && consistentGrid;
}

void Qstls::computeAdrIet() {
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  assert(adrOld.size() > 0);
  // Setup interpolators
  const Interpolator1D ssfi(wvg, ssfOld);
  const Interpolator1D bfi(wvg, bf);
  vector<Interpolator1D> dlfci(nl);
  Interpolator1D tmp(wvg, ssfOld);
  for (int l = 0; l < nl; ++l) {
    vector<double> dlfc(nx);
    for (int i = 0; i < nx; ++i) {
      dlfc[i] = (idr(i, l) > 0.0) ? adrOld(i, l) / idr(i, l) : 0;
    }
    dlfci[l].reset(wvg[0], dlfc[0], nx);
  }
  // Compute qstls-iet contribution to the adr
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  Vector2D adrIet(nx, nl);
  auto loopFunc = [&](int i) -> void {
    Integrator2D itgPrivate(in.getIntError());
    Vector3D adrFixedPrivate;
    readAdrFixedFile(adrFixedPrivate, adrFixedIetFileInfo.at(i).first, true);
    AdrIet adrTmp(in.getDegeneracy(),
                  wvg.front(),
                  wvg.back(),
                  wvg[i],
                  ssfi,
                  dlfci,
                  bfi,
                  itgGrid,
                  itgPrivate);
    adrTmp.get(wvg, adrFixedPrivate, adrIet);
  };
  const auto &loopData = parallelFor(loopFunc, nx, in.getNThreads());
  gatherLoopData(adrIet.data(), loopData, nl);
  // Sum qstls and qstls-iet contributions to adr
  adr.sum(adrIet);
}

void Qstls::computeAdrFixedIet() {
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  vector<int> idx;
  // Check which files have to be created
  getAdrFixedIetFileInfo();
  for (const auto &member : adrFixedIetFileInfo) {
    if (!member.second.second) { idx.push_back(member.first); };
  }
  // Check that all ranks found the same number of files
  int nFilesToWrite = idx.size();
  if (!isEqualOnAllRanks(idx.size())) {
    throwError("Not all ranks can access the files with the fixed "
               " component of the auxiliary density response");
  }
  if (nFilesToWrite == 0) { return; }
  // Barrier to ensure that all ranks have identified all the files that have to
  // be written
  barrier();
  // Write necessary files
  auto loopFunc = [&](int i) -> void {
    Integrator1D itgPrivate(in.getIntError());
    Vector3D res(nl, nx, nx);
    AdrFixedIet adrTmp(in.getDegeneracy(),
                       wvg.front(),
                       wvg.back(),
                       wvg[idx[i]],
                       mu,
                       itgPrivate);
    adrTmp.get(wvg, res);
    writeAdrFixedFile(res, adrFixedIetFileInfo.at(idx[i]).first);
  };
  parallelFor(loopFunc, nFilesToWrite, in.getNThreads());
  // Update content of adrFixedIetFileInfo
  getAdrFixedIetFileInfo();
}

void Qstls::getAdrFixedIetFileInfo() {
  const int nx = wvg.size();
  adrFixedIetFileInfo.clear();
  for (int i = 0; i < nx; ++i) {
    try {
      string name =
          fmt::format("adr_fixed_theta{:.3f}_matsubara{:}_{}_wv{:.5f}.bin",
                      in.getDegeneracy(),
                      in.getNMatsubara(),
                      in.getTheory(),
                      wvg[i]);
      if (!in.getFixedIet().empty()) {
        std::filesystem::path fullPath = in.getFixedIet();
        fullPath /= name;
        name = fullPath.string();
      }
      const bool found = std::filesystem::exists(name);
      const pair<string, bool> filePair = pair<string, bool>(name, found);
      adrFixedIetFileInfo.insert(pair<int, decltype(filePair)>(i, filePair));
    } catch (...) {
      throwError("Error in the output file for the fixed component"
                 " of the auxiliary density response.");
    }
  }
}

// Recovery files
void Qstls::writeRecovery() {
  if (!isRoot()) { return; }
  ofstream file;
  file.open(recoveryFileName, ios::binary);
  if (!file.is_open()) {
    throwError("Recovery file " + recoveryFileName + " could not be created.");
  }
  int nx = wvg.size();
  int nl = in.getNMatsubara();
  writeDataToBinary<int>(file, nx);
  writeDataToBinary<int>(file, nl);
  writeDataToBinary<double>(file, in.getDegeneracy());
  writeDataToBinary<vector<double>>(file, wvg);
  writeDataToBinary<vector<double>>(file, ssf);
  writeDataToBinary<Vector2D>(file, adr);
  writeDataToBinary<Vector3D>(file, adrFixed);
  file.close();
  if (!file) {
    throwError("Error in writing the recovery file " + recoveryFileName);
  }
}

void Qstls::readRecovery(vector<double> &wvg_,
                         vector<double> &ssf_,
                         Vector2D &adr_,
                         Vector3D &adrFixed_,
                         double &Theta,
                         int &nl) const {
  const string &fileName = in.getRecoveryFileName();
  if (fileName.empty()) { return; }
  ifstream file;
  file.open(fileName, ios::binary);
  if (!file.is_open()) {
    throwError("Input file " + fileName + " could not be opened.");
  }
  int nx;
  readDataFromBinary<int>(file, nx);
  readDataFromBinary<int>(file, nl);
  readDataFromBinary<double>(file, Theta);
  wvg_.resize(nx);
  ssf_.resize(nx);
  adr_.resize(nx, nl);
  adrFixed_.resize(nx, nl, nx);
  readDataFromBinary<vector<double>>(file, wvg_);
  readDataFromBinary<vector<double>>(file, ssf_);
  readDataFromBinary<Vector2D>(file, adr_);
  readDataFromBinary<Vector3D>(file, adrFixed_);
  file.close();
  if (!file) { throwError("Error in reading the file " + fileName); }
}

// -----------------------------------------------------------------
// AdrBase class
// -----------------------------------------------------------------

// Compute static structure factor
double AdrBase::ssf(const double &y) const { return ssfi.eval(y); }

// -----------------------------------------------------------------
// Adr class
// -----------------------------------------------------------------

// Compute fixed component
double Adr::fix(const double &y) const { return fixi.eval(y); }

// Integrand
double Adr::integrand(const double &y) const { return fix(y) * (ssf(y) - 1.0); }

// Get result of integration
void Adr::get(const vector<double> &wvg, const Vector3D &fixed, Vector2D &res) {
  const int nx = wvg.size();
  const int nl = fixed.size(1);
  auto it = lower_bound(wvg.begin(), wvg.end(), x);
  assert(it != wvg.end());
  size_t ix = distance(wvg.begin(), it);
  if (x == 0.0) {
    res.fill(ix, 0.0);
    return;
  }
  const auto itgParam = ItgParam(yMin, yMax);
  for (int l = 0; l < nl; ++l) {
    fixi.reset(wvg[0], fixed(ix, l), nx);
    auto func = [&](const double &y) -> double { return integrand(y); };
    itg.compute(func, itgParam);
    res(ix, l) = itg.getSolution();
    res(ix, l) *= (l == 0) ? isc0 : isc;
  }
}

// -----------------------------------------------------------------
// AdrFixed class
// -----------------------------------------------------------------

// Get fixed component
void AdrFixed::get(vector<double> &wvg, Vector3D &res) const {
  const int nx = wvg.size();
  const int nl = res.size(1);
  if (x == 0.0) { res.fill(0, 0.0); };
  const double x2 = x * x;
  auto it = find(wvg.begin(), wvg.end(), x);
  assert(it != wvg.end());
  size_t ix = distance(wvg.begin(), it);
  for (int l = 0; l < nl; ++l) {
    for (int i = 0; i < nx; ++i) {
      const double xq = x * wvg[i];
      auto tMin = x2 - xq;
      auto tMax = x2 + xq;
      auto func1 = [&](const double &q) -> double { return integrand1(q, l); };
      auto func2 = [&](const double &t) -> double {
        return integrand2(t, wvg[i], l);
      };
      itg.compute(func1, func2, Itg2DParam(qMin, qMax, tMin, tMax), itgGrid);
      res(ix, l, i) = itg.getSolution();
    }
  }
}

// Integrands for the fixed component
double AdrFixed::integrand1(const double &q, const double &l) const {
  if (l == 0)
    return q / (exp(q * q / Theta - mu) + exp(-q * q / Theta + mu) + 2.0);
  return q / (exp(q * q / Theta - mu) + 1.0);
}

double
AdrFixed::integrand2(const double &t, const double &y, const double &l) const {
  const double q = itg.getX();
  if (y == 0) { return 0.0; };
  const double x2 = x * x;
  const double y2 = y * y;
  const double q2 = q * q;
  const double txq = 2.0 * x * q;
  if (l == 0) {
    if (t == txq) { return 2.0 * q2 / (y2 + 2.0 * txq - x2); };
    if (x == y && t == 0.0) { return q; };
    const double t2 = t * t;
    double logarg = (t + txq) / (t - txq);
    logarg = (logarg < 0.0) ? -logarg : logarg;
    return y / (2.0 * t + y2 - x2)
           * ((q2 - t2 / (4.0 * x2)) * log(logarg) + q * t / x);
  }
  if (x == y && t == 0.0) { return 0.0; };
  const double tplT = 2.0 * M_PI * l * Theta;
  const double tplT2 = tplT * tplT;
  const double txqpt = txq + t;
  const double txqmt = txq - t;
  const double txqpt2 = txqpt * txqpt;
  const double txqmt2 = txqmt * txqmt;
  const double logarg = (txqpt2 + tplT2) / (txqmt2 + tplT2);
  return y / (2.0 * t + y * y - x * x) * log(logarg);
}

// -----------------------------------------------------------------
// AdrIet class
// -----------------------------------------------------------------

// Compute dynamic local field correction
double AdrIet::dlfc(const double &y, const int &l) const {
  return dlfci[l].eval(y);
}

// Compute auxiliary density response
double AdrIet::bf(const double &y) const { return bfi.eval(y); }

// Compute fixed component
double AdrIet::fix(const double &x, const double &y) const {
  return fixi.eval(x, y);
}

// Integrands
double AdrIet::integrand1(const double &q, const int &l) const {
  if (q == 0.0) { return 0.0; }
  const double p1 = (1 - bf(q)) * ssf(q);
  const double p2 = dlfc(q, l) * (ssf(q) - 1.0);
  return (p1 - p2 - 1.0) / q;
}

double AdrIet::integrand2(const double &y) const {
  const double q = itg.getX();
  return y * fix(q, y) * (ssf(y) - 1.0);
}

// Get result of integration
void AdrIet::get(const vector<double> &wvg,
                 const Vector3D &fixed,
                 Vector2D &res) {
  const int nx = wvg.size();
  const int nl = fixed.size(0);
  auto it = lower_bound(wvg.begin(), wvg.end(), x);
  assert(it != wvg.end());
  size_t ix = distance(wvg.begin(), it);
  if (x == 0.0) {
    res.fill(ix, 0.0);
    return;
  }
  for (int l = 0; l < nl; ++l) {
    fixi.reset(wvg[0], wvg[0], fixed(l), nx, nx);
    auto yMin = [&](const double &q) -> double {
      return (q > x) ? q - x : x - q;
    };
    auto yMax = [&](const double &q) -> double { return min(qMax, q + x); };
    auto func1 = [&](const double &q) -> double { return integrand1(q, l); };
    auto func2 = [&](const double &y) -> double { return integrand2(y); };
    itg.compute(func1, func2, Itg2DParam(qMin, qMax, yMin, yMax), itgGrid);
    res(ix, l) = itg.getSolution();
    res(ix, l) *= (l == 0) ? isc0 : isc;
  }
}

// -----------------------------------------------------------------
// AdrFixedIet class
// -----------------------------------------------------------------

// get fixed component
void AdrFixedIet::get(vector<double> &wvg, Vector3D &res) const {
  if (x == 0) {
    res.fill(0.0);
    return;
  }
  const int nx = wvg.size();
  const int nl = res.size(0);
  const auto itgParam = ItgParam(tMin, tMax);
  for (int l = 0; l < nl; ++l) {
    for (int i = 0; i < nx; ++i) {
      if (wvg[i] == 0.0) {
        res.fill(l, i, 0.0);
        continue;
      }
      for (int j = 0; j < nx; ++j) {
        auto func = [&](const double &t) -> double {
          return integrand(t, wvg[j], wvg[i], l);
        };
        itg.compute(func, itgParam);
        res(l, i, j) = itg.getSolution();
      }
    }
  }
}

// Integrand for the fixed component
double AdrFixedIet::integrand(const double &t,
                              const double &y,
                              const double &q,
                              const double &l) const {
  const double x2 = x * x;
  const double q2 = q * q;
  const double y2 = y * y;
  const double t2 = t * t;
  const double fxt = 4.0 * x * t;
  const double qmypx = q2 - y2 + x2;
  if (l == 0) {
    double logarg = (qmypx + fxt) / (qmypx - fxt);
    if (logarg < 0.0) logarg = -logarg;
    return t / (exp(t2 / Theta - mu) + exp(-t2 / Theta + mu) + 2.0)
           * ((t2 - qmypx * qmypx / (16.0 * x2)) * log(logarg)
              + (t / x) * qmypx / 2.0);
  }
  const double fplT = 4.0 * M_PI * l * Theta;
  const double fplT2 = fplT * fplT;
  const double logarg = ((qmypx + fxt) * (qmypx + fxt) + fplT2)
                        / ((qmypx - fxt) * (qmypx - fxt) + fplT2);
  return t / (exp(t2 / Theta - mu) + 1.0) * log(logarg);
}

// -----------------------------------------------------------------
// AdrGround class
// -----------------------------------------------------------------

// Get
double AdrGround::get() {
  auto tMin = [&](const double &y) -> double { return x * (x + y); };
  auto tMax = [&](const double &y) -> double { return x * (x - y); };
  auto func1 = [&](const double &y) -> double { return integrand1(y); };
  auto func2 = [&](const double &t) -> double { return integrand2(t); };
  itg.compute(func1, func2, Itg2DParam(0, yMax, tMin, tMax), {});
  return -(3.0 / 8.0) * itg.getSolution();
}

double AdrGround::integrand1(const double &y) const {
  return y * (ssfi.eval(y) - 1.0);
}

double AdrGround::integrand2(const double &t) const {
  if (x == 0.0) { return 0.0; }
  const double y = itg.getX();
  const double x2 = x * x;
  const double Omega2 = Omega * Omega;
  const double t2 = t * t;
  const double y2 = y * y;
  const double tx = 2.0 * x;
  const double tptx = t + tx;
  const double tmtx = t - tx;
  const double tptx2 = tptx * tptx;
  const double tmtx2 = tmtx * tmtx;
  const double logarg = (tptx2 + Omega2) / (tmtx2 + Omega2);
  const double part1 =
      (0.5 - t2 / (8.0 * x2) + Omega2 / (8.0 * x2)) * log(logarg);
  const double part2 =
      0.5 * Omega * t / x2 * (atan(tptx / Omega) - atan(tmtx / Omega));
  const double part3 = part1 - part2 + t / x;
  return part3 / (2.0 * t + y2 - x2);
}

// -----------------------------------------------------------------
// Qssf class
// -----------------------------------------------------------------

double Qssf::get() const {
  if (rs == 0.0) return ssfHF;
  if (x == 0.0) return 0.0;
  const double f2 = 1 - bf;
  double f3 = 0.0;
  for (int l = 0; l < nl; ++l) {
    const double f4 = f2 * idr[l];
    const double f5 = 1.0 + ip * (f4 - adr[l]);
    const double f6 = idr[l] * (f4 - adr[l]) / f5;
    f3 += (l == 0) ? f6 : 2 * f6;
  }
  return ssfHF - 1.5 * ip * Theta * f3;
}

// -----------------------------------------------------------------
// QssfGround class
// -----------------------------------------------------------------

double QssfGround::get() {
  if (x == 0.0) return 0.0;
  if (rs == 0.0) return ssfHF;
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg.compute(func, ItgParam(0, OmegaMax));
  return 1.5 / (M_PI)*itg.getSolution() + ssfHF;
}

double QssfGround::integrand(const double &Omega) const {
  Integrator2D itg2 = Integrator2D(itg.getAccuracy());
  const double ip = 4.0 * lambda * rs / (M_PI * x * x);
  const double idr = IdrGround(x, Omega).get();
  const double adr = AdrGround(x, Omega, ssfi, xMax, itg2).get();
  return idr / (1.0 + ip * (idr - adr)) - idr;
}
