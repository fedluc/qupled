#include "qstls.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "vector_util.hpp"
#include <SQLiteCpp/SQLiteCpp.h>

using namespace std;
using namespace vecUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;

// -----------------------------------------------------------------
// QSTLS class
// -----------------------------------------------------------------

Qstls::Qstls(const QstlsInput &in_, const bool verbose_)
    : Stls(in_, verbose_),
      in(in_) {
  // Set name for the fixed adr output file
  adrFixedDatabaseName = fmt::format("{}", in.getTheory());
  // Allocate arrays
  const size_t nx = wvg.size();
  const size_t nl = in.getNMatsubara();
  adr.resize(nx, nl);
  ssfNew.resize(nx);
  ssfOld.resize(nx);
  adrFixed.resize(nx, nl, nx);
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
}

// qstls iterations
void Qstls::doIterations() {
  const int maxIter = in.getNIter();
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
  // From guess in input
  if (initialGuessFromInput()) { return; }
  // Default
  Rpa rpa(in, false);
  int status = rpa.compute();
  if (status != 0) {
    throwError("Failed to compute the default initial guess");
  }
  ssfOld = rpa.getSsf();
}

bool Qstls::initialGuessFromInput() {
  const auto &guess = in.getGuess();
  const int nx = wvg.size();
  const Interpolator1D ssfi(guess.wvg, guess.ssf);
  if (!ssfi.isValid()) { return false; }
  const double xMax = guess.wvg.back();
  for (int i = 0; i < nx; ++i) {
    const double &x = wvg[i];
    ssfOld[i] = (x <= xMax) ? ssfi.eval(x) : 1.0;
  }
  return true;
}

// Compute auxiliary density response
void Qstls::computeAdr() {
  if (in.getDegeneracy() == 0.0) { return; }
  const int nx = wvg.size();
  const shared_ptr<Interpolator1D> ssfi =
      make_shared<Interpolator1D>(wvg, ssfOld);
  for (int i = 0; i < nx; ++i) {
    QstlsUtil::Adr adrTmp(
        in.getDegeneracy(), wvg.front(), wvg.back(), wvg[i], ssfi, itg);
    adrTmp.get(wvg, adrFixed, adr);
  }
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
  }
  const double Theta = in.getDegeneracy();
  const double rs = in.getCoupling();
  const int nx = wvg.size();
  for (int i = 0; i < nx; ++i) {
    QstlsUtil::Ssf ssfTmp(wvg[i], Theta, rs, ssfHF[i], idr[i], adr[i]);
    ssfNew[i] = ssfTmp.get();
  }
}

// Compute static structure factor at zero temperature
void Qstls::computeSsfGround() {
  const shared_ptr<Interpolator1D> ssfi =
      make_shared<Interpolator1D>(wvg, ssfOld);
  const double rs = in.getCoupling();
  const double OmegaMax = in.getFrequencyCutoff();
  const size_t nx = wvg.size();
  const double xMax = wvg.back();
  auto loopFunc = [&](int i) -> void {
    shared_ptr<Integrator1D> itgTmp = make_shared<Integrator1D>(*itg);
    QstlsUtil::SsfGround ssfTmp(
        wvg[i], rs, ssfHF[i], xMax, OmegaMax, ssfi, itgTmp);
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
}

void Qstls::computeAdrFixed() {
  if (in.getDegeneracy() == 0.0) { return; }
  // Check if it adrFixed can be loaded from input
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  if (in.getFixedRunId() != DEFAULT_INT) {
    adrFixed.resize(nx, nl, nx);
    readAdrFixed(adrFixed, adrFixedDatabaseName, in.getFixedRunId());
    return;
  }
  // Compute from scratch
  fflush(stdout);
  const int nxnl = nx * nl;
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  // Parallel for loop (Hybrid MPI and OpenMP)
  auto loopFunc = [&](int i) -> void {
    shared_ptr<Integrator2D> itg2 = make_shared<Integrator2D>(in.getIntError());
    QstlsUtil::AdrFixed adrTmp(
        in.getDegeneracy(), wvg.front(), wvg.back(), wvg[i], mu, itgGrid, itg2);
    adrTmp.get(wvg, adrFixed);
  };
  const auto &loopData = parallelFor(loopFunc, nx, in.getNThreads());
  gatherLoopData(adrFixed.data(), loopData, nxnl);
  // Write result to output file
  if (isRoot()) { writeAdrFixed(adrFixed, adrFixedDatabaseName); }
}

void Qstls::writeAdrFixed(const Vector3D &res, const string &name) const {
  try {
    DatabaseInfo dbInfo = in.getDatabaseInfo();
    SQLite::Database db(dbInfo.name,
                        SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE);
    // Create table if it doesn't exist
    const string createTable = fmt::format(QstlsUtil::SQL_CREATE_TABLE,
                                           QstlsUtil::SQL_TABLE_NAME,
                                           dbInfo.runTableName);
    db.exec(createTable);
    // Prepare binary data
    const void *adrData = static_cast<const void *>(res.data());
    const int adrSize = static_cast<int>(res.size() * sizeof(double));
    // Insert data
    const string insert =
        fmt::format(QstlsUtil::SQL_INSERT, QstlsUtil::SQL_TABLE_NAME);
    SQLite::Statement statement(db, insert);
    statement.bind(1, dbInfo.runId);
    statement.bind(2, name);
    statement.bind(3, adrData, adrSize);
    statement.exec();
  } catch (const std::exception &e) {
    throwError("Failed to write to database: " + string(e.what()));
  }
}

void Qstls::readAdrFixed(Vector3D &res, const string &name, int runId) const {
  try {
    DatabaseInfo dbInfo = in.getDatabaseInfo();
    SQLite::Database db(dbInfo.name, SQLite::OPEN_READONLY);
    const string select =
        fmt::format(QstlsUtil::SQL_SELECT, QstlsUtil::SQL_TABLE_NAME);
    SQLite::Statement statement(db, select);
    statement.bind(1, runId);
    statement.bind(2, name);
    if (statement.executeStep()) {
      const void *adrData = statement.getColumn(0).getBlob();
      int adrBytes = statement.getColumn(0).getBytes();
      if (static_cast<size_t>(adrBytes) != res.size() * sizeof(double)) {
        throwError(fmt::format("Size mismatch: expected {} bytes, got {} bytes",
                               res.size() * sizeof(double),
                               adrBytes));
      }
      std::memcpy(res.data(), adrData, adrBytes);
    }
  } catch (const std::exception &e) {
    throwError("Failed to read from database: " + string(e.what()));
  }
}

// -----------------------------------------------------------------
// AdrBase class
// -----------------------------------------------------------------

// Compute static structure factor
double QstlsUtil::AdrBase::ssf(const double &y) const { return ssfi->eval(y); }

// -----------------------------------------------------------------
// Adr class
// -----------------------------------------------------------------

// Compute fixed component
double QstlsUtil::Adr::fix(const double &y) const { return fixi.eval(y); }

// Integrand
double QstlsUtil::Adr::integrand(const double &y) const {
  return fix(y) * (ssf(y) - 1.0);
}

// Get result of integration
void QstlsUtil::Adr::get(const vector<double> &wvg,
                         const Vector3D &fixed,
                         Vector2D &res) {
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
    fixi.reset(wvg[0], fixed(ix, l, 0), nx);
    auto func = [&](const double &y) -> double { return integrand(y); };
    itg->compute(func, itgParam);
    res(ix, l) = itg->getSolution();
    res(ix, l) *= (l == 0) ? isc0 : isc;
  }
}

// -----------------------------------------------------------------
// AdrFixed class
// -----------------------------------------------------------------

// Get fixed component
void QstlsUtil::AdrFixed::get(const vector<double> &wvg, Vector3D &res) const {
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
      itg->compute(func1, func2, Itg2DParam(qMin, qMax, tMin, tMax), itgGrid);
      res(ix, l, i) = itg->getSolution();
    }
  }
}

// Integrands for the fixed component
double QstlsUtil::AdrFixed::integrand1(const double &q, const double &l) const {
  if (l == 0)
    return q / (exp(q * q / Theta - mu) + exp(-q * q / Theta + mu) + 2.0);
  return q / (exp(q * q / Theta - mu) + 1.0);
}

double QstlsUtil::AdrFixed::integrand2(const double &t,
                                       const double &y,
                                       const double &l) const {
  const double q = itg->getX();
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
// AdrGround class
// -----------------------------------------------------------------

// Get
double QstlsUtil::AdrGround::get() {
  auto tMin = [&](const double &y) -> double { return x * (x + y); };
  auto tMax = [&](const double &y) -> double { return x * (x - y); };
  auto func1 = [&](const double &y) -> double { return integrand1(y); };
  auto func2 = [&](const double &t) -> double { return integrand2(t); };
  itg->compute(func1, func2, Itg2DParam(0, yMax, tMin, tMax), {});
  return -(3.0 / 8.0) * itg->getSolution();
}

double QstlsUtil::AdrGround::integrand1(const double &y) const {
  return y * (ssfi->eval(y) - 1.0);
}

double QstlsUtil::AdrGround::integrand2(const double &t) const {
  if (x == 0.0) { return 0.0; }
  const double y = itg->getX();
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

double QstlsUtil::Ssf::get() const {
  if (rs == 0.0) return ssfHF;
  if (x == 0.0) return 0.0;
  double suml = 0.0;
  for (size_t l = 0; l < idr.size(); ++l) {
    const double &idrl = idr[l];
    const double &adrl = adr[l];
    const double denom = 1.0 + ip * (idrl - adrl);
    const double f = idrl * (idrl - adrl) / denom;
    suml += (l == 0) ? f : 2 * f;
  }
  return ssfHF - 1.5 * ip * Theta * suml;
}

// -----------------------------------------------------------------
// QssfGround class
// -----------------------------------------------------------------

double QstlsUtil::SsfGround::get() {
  if (x == 0.0) return 0.0;
  if (rs == 0.0) return ssfHF;
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg->compute(func, ItgParam(0, OmegaMax));
  return 1.5 / (M_PI)*itg->getSolution() + ssfHF;
}

double QstlsUtil::SsfGround::integrand(const double &Omega) const {
  shared_ptr<Integrator2D> itg2 = make_shared<Integrator2D>(itg->getAccuracy());
  const double idr = HFUtil::IdrGround(x, Omega).get();
  const double adr = QstlsUtil::AdrGround(x, Omega, ssfi, xMax, itg2).get();
  return idr / (1.0 + ip * (idr - adr)) - idr;
}
