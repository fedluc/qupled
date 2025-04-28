#include "qstlsiet.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "stlsiet.hpp"
#include "vector_util.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <filesystem>
#include <fmt/core.h>
#include <numeric>

using namespace std;
using namespace vecUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;
using ItgType = Integrator1D::Type;

// -----------------------------------------------------------------
// QSTLS-IET class
// -----------------------------------------------------------------

QstlsIet::QstlsIet(const QstlsIetInput &in_)
    : Qstls(in_, true),
      Iet(in_, in_, true),
      in(in_) {
  // Throw error message for ground state calculations
  if (in.getDegeneracy() == 0.0) {
    throwError("Ground state calculations are not available "
               "for the quantum IET schemes");
  }
  // Allocate arrays
  adrOld.resize(adr.size(0), adr.size(1));
}

int QstlsIet::compute() {
  try {
    init();
    Qstls::println("Structural properties calculation ...");
    doIterations();
    Qstls::println("Done");
    return 0;
  } catch (const runtime_error &err) {
    cerr << err.what() << endl;
    return 1;
  }
}

void QstlsIet::init() {
  Qstls::init();
  Iet::init();
  Qstls::print(
      "Computing fixed component of the iet auxiliary density response: ");
  computeAdrFixed();
  Qstls::println("Done");
}

void QstlsIet::doIterations() {
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
    Qstls::println(fmt::format("--- iteration {:d} ---", counter));
    Qstls::println(fmt::format("Elapsed time: {:.3f} seconds", toc - tic));
    Qstls::println(fmt::format("Residual error: {:.5e}", err));
    fflush(stdout);
  }
  // Set static structure factor for output
  ssf = ssfOld;
}

void QstlsIet::initialGuess() {
  Qstls::initialGuess();
  assert(!adrOld.empty());
  if (initialGuessFromInput()) { return; }
  adrOld.fill(0.0);
}

bool QstlsIet::initialGuessFromInput() {
  const auto &guess = in.getGuess();
  const vector<double> &wvg = in.getWaveVectorGrid();
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const int nx_ = guess.adr.size(0);
  const int nl_ = guess.adr.size(1);
  const double xMax = (guess.wvg.empty()) ? 0.0 : guess.wvg.back();
  vector<Interpolator1D> itp(nl_);
  for (int l = 0; l < nl_; ++l) {
    vector<double> tmp(nx_);
    for (int i = 0; i < nx_; ++i) {
      tmp[i] = guess.adr(i, l);
    }
    itp[l].reset(guess.wvg[0], tmp[0], nx_);
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

void QstlsIet::computeSsf() {
  const double Theta = in.getDegeneracy();
  const double rs = in.getCoupling();
  const vector<double> &wvg = in.getWaveVectorGrid();
  const int nx = wvg.size();
  const int nl = idr.size(1);
  const vector<double> &bf = getBf();
  for (int i = 0; i < nx; ++i) {
    const double bfi = bf[i];
    QstlsIetUtil::Ssf ssfTmp(
        wvg[i], Theta, rs, ssfHF[i], nl, &idr(i), &adr(i), bfi);
    ssfNew[i] = ssfTmp.get();
  }
}

void QstlsIet::updateSolution() {
  Qstls::updateSolution();
  const double aMix = in.getMixingParameter();
  adrOld.mult(1 - aMix);
  adrOld.linearCombination(adr, aMix);
}

void QstlsIet::computeAdr() {
  Qstls::computeAdr();
  const vector<double> &wvg = in.getWaveVectorGrid();
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  assert(adrOld.size() > 0);
  // Setup interpolators
  const shared_ptr<Interpolator1D> ssfi =
      make_shared<Interpolator1D>(wvg, ssfOld);
  const shared_ptr<Interpolator1D> bfi =
      make_shared<Interpolator1D>(wvg, getBf());
  vector<shared_ptr<Interpolator1D>> dlfci(nl);
  shared_ptr<Interpolator1D> tmp = make_shared<Interpolator1D>(wvg, ssfOld);
  for (int l = 0; l < nl; ++l) {
    vector<double> dlfc(nx);
    for (int i = 0; i < nx; ++i) {
      dlfc[i] = (idr(i, l) > 0.0) ? adrOld(i, l) / idr(i, l) : 0;
    }
    dlfci[l] = std::make_shared<Interpolator1D>(wvg, dlfc);
  }
  // Compute qstls-iet contribution to the adr
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  Vector2D adrIet(nx, nl);
  auto loopFunc = [&](int i) -> void {
    shared_ptr<Integrator2D> itgPrivate =
        make_shared<Integrator2D>(in.getIntError());
    Vector3D adrFixedPrivate(nl, nx, nx);
    const string name = fmt::format("{}_{:d}", in.getTheory(), i);
    const int runId = (in.getFixedRunId() != DEFAULT_INT)
                          ? in.getFixedRunId()
                          : in.getDatabaseInfo().runId;
    readAdrFixed(adrFixedPrivate, name, runId);
    QstlsIetUtil::AdrIet adrTmp(in.getDegeneracy(),
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
  // Compute static local field correction
  for (int i = 0; i < nx; ++i) {
    slfc[i] = adr(i, 0) / idr(i, 0);
  };
}

void QstlsIet::computeAdrFixed() {
  if (in.getFixedRunId() != DEFAULT_INT) { return; }
  const vector<double> &wvg = in.getWaveVectorGrid();
  const int nx = wvg.size();
  const int nl = in.getNMatsubara();
  const double &xStart = wvg.front();
  const double &xEnd = wvg.back();
  const double &theta = in.getDegeneracy();
  for (const auto &x : wvg) {
    Vector3D res(nl, nx, nx);
    auto loopFunc = [&](int l) -> void {
      auto itgPrivate = make_shared<Integrator1D>(in.getIntError());
      QstlsIetUtil::AdrFixedIet adrTmp(theta, xStart, xEnd, x, mu, itgPrivate);
      adrTmp.get(l, wvg, res);
    };
    const auto &loopData = parallelFor(loopFunc, nl, in.getNThreads());
    gatherLoopData(res.data(), loopData, nx * nx);
    if (isRoot()) {
      const size_t idx = distance(wvg.begin(), find(wvg.begin(), wvg.end(), x));
      const string name = fmt::format("{}_{:d}", in.getTheory(), idx);
      writeAdrFixed(res, name);
    }
  }
  // TODO: Check that all ranks can access the database
  barrier();
}

// -----------------------------------------------------------------
// AdrIet class
// -----------------------------------------------------------------

// Compute dynamic local field correction
double QstlsIetUtil::AdrIet::dlfc(const double &y, const int &l) const {
  return dlfci[l]->eval(y);
}

// Compute auxiliary density response
double QstlsIetUtil::AdrIet::bf(const double &y) const { return bfi->eval(y); }

// Compute fixed component
double QstlsIetUtil::AdrIet::fix(const double &x, const double &y) const {
  return fixi.eval(x, y);
}

// Integrands
double QstlsIetUtil::AdrIet::integrand1(const double &q, const int &l) const {
  if (q == 0.0) { return 0.0; }
  const double p1 = (1 - bf(q)) * ssf(q);
  const double p2 = dlfc(q, l) * (ssf(q) - 1.0);
  return (p1 - p2 - 1.0) / q;
}

double QstlsIetUtil::AdrIet::integrand2(const double &y) const {
  const double q = itg->getX();
  return y * fix(q, y) * (ssf(y) - 1.0);
}

// Get result of integration
void QstlsIetUtil::AdrIet::get(const vector<double> &wvg,
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
    itg->compute(func1, func2, Itg2DParam(qMin, qMax, yMin, yMax), itgGrid);
    res(ix, l) = itg->getSolution();
    res(ix, l) *= (l == 0) ? isc0 : isc;
  }
}

// -----------------------------------------------------------------
// AdrFixedIet class
// -----------------------------------------------------------------

// get fixed component
void QstlsIetUtil::AdrFixedIet::get(int l,
                                    const vector<double> &wvg,
                                    Vector3D &res) const {
  if (x == 0.0) {
    res.fill(l, 0.0);
    return;
  }
  const int nx = wvg.size();
  const auto itgParam = ItgParam(tMin, tMax);
  for (int i = 0; i < nx; ++i) {
    if (wvg[i] == 0.0) {
      res.fill(l, i, 0.0);
      continue;
    }
    for (int j = 0; j < nx; ++j) {
      auto func = [&](const double &t) -> double {
        return integrand(t, wvg[j], wvg[i], l);
      };
      itg->compute(func, itgParam);
      res(l, i, j) = itg->getSolution();
    }
  }
}

// Integrand for the fixed component
double QstlsIetUtil::AdrFixedIet::integrand(const double &t,
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
// Ssf class
// -----------------------------------------------------------------

double QstlsIetUtil::Ssf::get() const {
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
