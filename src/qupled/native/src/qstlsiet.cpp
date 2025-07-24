#include "qstlsiet.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "stlsiet.hpp"
#include "vector_util.hpp"
#include <SQLiteCpp/SQLiteCpp.h>
#include <filesystem>
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

QstlsIet::QstlsIet(const std::shared_ptr<const QstlsIetInput> &in_)
    : Qstls(in_, true),
      iet(in_, wvg) {
  // Throw error message for ground state calculations
  if (in().getDegeneracy() == 0.0) {
    throwError(
        "Ground state calculations are not implemented for this scheme.");
  }
  // Allocate grid for the 2D integrator
  const bool segregatedItg = in().getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  // Allocate arrays
  const size_t nx = wvg.size();
  const size_t nl = in().getNMatsubara();
  lfcIet.resize(nx, nl);
}

void QstlsIet::init() {
  Qstls::init();
  print("Computing fixed component of the iet auxiliary density response: ");
  computeAdrFixed();
  println("Done");
  iet.init();
}

bool QstlsIet::initialGuessFromInput() {
  const bool ssfIsSetFromInput = Qstls::initialGuessFromInput();
  if (!ssfIsSetFromInput) { return false; }
  const bool lfcIsSetFromInput = iet.initialGuessFromInput(lfc);
  if (!lfcIsSetFromInput) { return false; }
  return true;
}

void QstlsIet::computeLfc() {
  const int nx = wvg.size();
  const int nl = in().getNMatsubara();
  // Setup interpolators
  const shared_ptr<Interpolator1D> ssfi = make_shared<Interpolator1D>(wvg, ssf);
  const shared_ptr<Interpolator1D> bfi =
      make_shared<Interpolator1D>(wvg, getBf());
  vector<shared_ptr<Interpolator1D>> lfci(nl);
  for (int l = 0; l < nl; ++l) {
    vector<double> lfcColumn(nx);
    for (int i = 0; i < nx; ++i) {
      lfcColumn[i] = lfc(i, l);
    }
    lfci[l] = std::make_shared<Interpolator1D>(wvg, lfcColumn);
  }
  // Compute the qstls contribution to the adr
  Qstls::computeLfc();
  // Compute qstls-iet contribution to the adr
  auto loopFunc = [&](int i) -> void {
    shared_ptr<Integrator2D> itgPrivate =
        make_shared<Integrator2D>(in().getIntError());
    Vector3D adrFixedPrivate(nl, nx, nx);
    const string name = formatUtil::format("{}_{:d}", in().getTheory(), i);
    const int runId = (in().getFixedRunId() != DEFAULT_INT)
                          ? in().getFixedRunId()
                          : in().getDatabaseInfo().runId;
    readAdrFixed(adrFixedPrivate, name, runId);
    QstlsIetUtil::AdrIet adrTmp(in().getDegeneracy(),
                                wvg.front(),
                                wvg.back(),
                                wvg[i],
                                ssfi,
                                lfci,
                                bfi,
                                itgGrid,
                                itgPrivate);
    adrTmp.get(wvg, adrFixedPrivate, lfcIet);
  };
  const auto &loopData = parallelFor(loopFunc, nx, in().getNThreads());
  gatherLoopData(lfcIet.data(), loopData, nl);
  // Sum qstls and qstls-iet contributions to the local field correction
  lfcIet.div(idr);
  lfcIet.fill(0, 0.0);
  lfc.sum(lfcIet);
  // Add the bridge function contribution
  const vector<double> &bf = getBf();
  for (int i = 0; i < nx; ++i) {
    span<double> lfcRow = lfc[i];
    std::for_each(
        lfcRow.begin(), lfcRow.end(), [&](double &vi) { vi += bf[i]; });
  }
}

void QstlsIet::computeAdrFixed() {
  if (in().getFixedRunId() != DEFAULT_INT) { return; }
  const int nx = wvg.size();
  const int nl = in().getNMatsubara();
  const double &xStart = wvg.front();
  const double &xEnd = wvg.back();
  const double &theta = in().getDegeneracy();
  for (const auto &x : wvg) {
    Vector3D res(nl, nx, nx);
    auto loopFunc = [&](int l) -> void {
      auto itgPrivate = make_shared<Integrator1D>(in().getIntError());
      QstlsIetUtil::AdrFixedIet adrTmp(theta, xStart, xEnd, x, mu, itgPrivate);
      adrTmp.get(l, wvg, res);
    };
    const auto &loopData = parallelFor(loopFunc, nl, in().getNThreads());
    gatherLoopData(res.data(), loopData, nx * nx);
    if (isRoot()) {
      const size_t idx = distance(wvg.begin(), find(wvg.begin(), wvg.end(), x));
      const string name = formatUtil::format("{}_{:d}", in().getTheory(), idx);
      writeAdrFixed(res, name);
    }
  }
  // Check that all ranks can access the database
  barrier();
  DatabaseInfo dbInfo = in().getDatabaseInfo();
  bool dbExists = std::filesystem::exists(dbInfo.name);
  if (!dbExists) {
    throwError(formatUtil::format(
        "Not all ranks can access the database file {}", dbInfo.name));
  }
}

// -----------------------------------------------------------------
// AdrIet class
// -----------------------------------------------------------------

// Compute dynamic local field correction
double QstlsIetUtil::AdrIet::lfc(const double &y, const int &l) const {
  return lfci[l]->eval(y);
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
  const double lfcql = lfc(q, l);
  const double ssfq = ssf(q);
  const double p1 = (1 - lfcql) * ssfq;
  const double p2 = -lfcql + bf(q);
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
    fixi.reset(wvg[0], wvg[0], fixed(l, 0, 0), nx, nx);
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