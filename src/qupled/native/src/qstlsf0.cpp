#include "qstlsf0.hpp"
#include "esa.hpp"
#include "mpi_util.hpp"
#include "num_util.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>

using namespace std;
using namespace MPIUtil;
using Itg2DParam = Integrator2D::Param;

namespace {

  constexpr double TOL = 1.0e-8;
  constexpr double DEFAULT_A_CUTOFF = 0.5;

} // namespace

// -----------------------------------------------------------------
// QSTLS-F0 class
// -----------------------------------------------------------------

QstlsF0::QstlsF0(const std::shared_ptr<const QstlsPimcInput> &in_)
    : Qstls(in_, true) {
  if (in().getDegeneracy() <= 0.0) {
    throwError("QSTLS-F0 is implemented only for finite degeneracy (theta > 0)");
  }
}

void QstlsF0::init() {
  Stls::init();
  print("Computing fixed component of the auxiliary density response (F0): ");
  fflush(stdout);
  computeAdrFixedF0();
  println("Done");
}

void QstlsF0::computeAdrFixedF0() {
  if (in().getDegeneracy() == 0.0) { return; }
  const int nx = wvg.size();
  const int nl = in().getNMatsubara();
  const QstlsF0Util::Distribution f0(in().getCoupling(),
                                     in().getDegeneracy(),
                                     mu,
                                     in().getPimcEta(),
                                     in().getPimcYSec(),
                                     in().getPimcACutoff());
  f0Grid = wvg;
  f0Values.resize(wvg.size());
  for (size_t i = 0; i < wvg.size(); ++i) { f0Values[i] = f0.value(wvg[i]); }

  if (in().getFixedRunId() != DEFAULT_INT) {
    adrFixed.resize(nx, nl, nx);
    readAdrFixed(adrFixed, adrFixedDatabaseName, in().getFixedRunId());
    return;
  }

  const int nxnl = nx * nl;
  const vector<double> itgGrid = wvg;
  atomic<int> completed{0};
  const int progressStep = max(1, nx / 10);

  auto loopFunc = [&](int i) -> void {
    shared_ptr<Integrator2D> itg2 = make_shared<Integrator2D>(in().getIntError());
    QstlsF0Util::AdrFixedF0 adrTmp(in().getDegeneracy(),
                                   wvg.front(),
                                   wvg.back(),
                                   wvg[i],
                                   mu,
                                   f0,
                                   itgGrid,
                                   itg2);
    adrTmp.get(wvg, adrFixed);
    const int done = ++completed;
    if (isRoot() && (done % progressStep == 0 || done == nx)) {
      cerr << "[QSTLS-F0] fixed ADR progress: " << done << "/" << nx << endl;
    }
  };

  const auto &loopData = parallelFor(loopFunc, nx, in().getNThreads());
  gatherLoopData(adrFixed.data(), loopData, nxnl);
  if (isRoot()) { writeAdrFixed(adrFixed, adrFixedDatabaseName); }
}

void QstlsF0::computeSsfFinite() {
  const double theta = in().getDegeneracy();
  const double rs = in().getCoupling();
  const size_t nx = wvg.size();
  const size_t nl = idr.size(1);
  for (size_t i = 0; i < nx; ++i) {
    const double x = wvg[i];
    if (x == 0.0) {
      ssf[i] = 0.0;
      continue;
    }
    const double ip = 4.0 * numUtil::lambda * rs / (M_PI * x * x);
    double suml = 0.0;
    for (size_t l = 0; l < nl; ++l) {
      const double phi = idr(i, l);
      const double g = lfc(i, l);
      const double denom = 1.0 + ip * phi * (1.0 - g);
      const double weight = (l == 0) ? 1.0 : 2.0;
      suml += weight * phi / denom;
    }
    ssf[i] = 1.5 * theta * suml;
  }
}

// -----------------------------------------------------------------
// Distribution helper
// -----------------------------------------------------------------

QstlsF0Util::Distribution::Distribution(const double &rs_,
                                        const double &theta_,
                                        const double &mu_,
                                        const double &etaIn,
                                        const double &ySecIn,
                                        const double &aCutoffIn)
    : rs(rs_),
      theta(theta_),
      mu(mu_),
      eta(0.0),
      ySec(0.0),
      aCutoff(0.0),
      tailCoeff(0.0) {
  const bool hasEta = !std::isnan(etaIn);
  const bool hasYSec = !std::isnan(ySecIn);
  const bool hasACutoff = !std::isnan(aCutoffIn);
  eta = hasEta ? etaIn : 8.0;
  aCutoff = hasACutoff ? aCutoffIn : DEFAULT_A_CUTOFF;

  const double g0 = ESAUtil::Slfc(rs, theta).onTop();
  const double rs2 = rs * rs;
  const double lambda2 = numUtil::lambda * numUtil::lambda;
  tailCoeff = 8.0 / (9.0 * M_PI * M_PI) * lambda2 * rs2 * std::abs(g0);

  ySec = hasYSec ? ySecIn : determineYSec();
  if (isRoot()) {
    cerr << "[QSTLS-F0] normalization moment: " << moment(eta, ySec)
         << " target=" << (1.0 / 3.0) << endl;
  }
}

double QstlsF0Util::Distribution::fermiDiracValue(const double &y) const {
  const double y2 = y * y;
  return 1.0 / (exp(y2 / theta - mu) + 1.0);
}

double QstlsF0Util::Distribution::fermiDiracDerivative(const double &y) const {
  const double f = fermiDiracValue(y);
  return -(2.0 * y / theta) * f * (1.0 - f);
}

double QstlsF0Util::Distribution::asymptoticValue(const double &y) const {
  const double yc = max(y, 1.0e-8);
  const double y2 = yc * yc;
  const double y4 = y2 * y2;
  return tailCoeff / (y4 * y4);
}

double QstlsF0Util::Distribution::asymptoticDerivative(const double &y) const {
  const double yc = max(y, 1.0e-8);
  return -8.0 * tailCoeff / pow(yc, 9);
}

double QstlsF0Util::Distribution::switchFunction(const double &y) const {
  if (y <= aCutoff) { return 0.0; }
  return 0.5 * (1.0 + tanh(eta * (y - ySec)));
}

double QstlsF0Util::Distribution::switchDerivative(const double &y) const {
  if (y <= aCutoff) { return 0.0; }
  const double arg = eta * (y - ySec);
  const double c = cosh(arg);
  return 0.5 * eta / (c * c);
}

double QstlsF0Util::Distribution::value(const double &y) const {
  return value(y, eta, ySec);
}

double QstlsF0Util::Distribution::value(const double &y,
                                        const double &etaLocal,
                                        const double &ySecLocal) const {
  const double ff = fermiDiracValue(y);
  const double fa = asymptoticValue(y);
  if (y <= aCutoff) { return ff; }
  const double a = 0.5 * (1.0 + tanh(etaLocal * (y - ySecLocal)));
  return (1.0 - a) * ff + a * fa;
}

double QstlsF0Util::Distribution::determineYSec() const { return solveYSec(eta); }

double QstlsF0Util::Distribution::derivative(const double &y) const {
  const double ff = fermiDiracValue(y);
  const double fa = asymptoticValue(y);
  const double dff = fermiDiracDerivative(y);
  const double dfa = asymptoticDerivative(y);
  const double a = switchFunction(y);
  const double da = switchDerivative(y);
  return (1.0 - a) * dff + a * dfa + da * (fa - ff);
}

double QstlsF0Util::Distribution::moment(const double &etaLocal,
                                         const double &ySecLocal) const {
  const double yUpper = 100.0;
  const int n = 8000;
  const double h = yUpper / static_cast<double>(n);
  double sum = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double y = i * h;
    const double w = (i == 0 || i == n) ? 0.5 : 1.0;
    sum += w * y * y * value(y, etaLocal, ySecLocal);
  }
  const double tail = tailCoeff / (5.0 * pow(yUpper, 5));
  return h * sum + tail;
}

double QstlsF0Util::Distribution::solveYSec(const double &etaLocal) const {
  const double target = 1.0 / 3.0;
  auto f = [&](double ySecLocal) { return moment(etaLocal, ySecLocal) - target; };
  double a = 0.01;
  double b = 50.0;
  double fa = f(a);
  double fb = f(b);
  double bestY = a;
  double bestErr = abs(fa);
  if (fa * fb > 0.0) {
    const int nScan = 240;
    double prevX = a;
    double prevF = fa;
    for (int i = 1; i <= nScan; ++i) {
      const double x = a + (b - a) * static_cast<double>(i) / nScan;
      const double fx = f(x);
      if (abs(fx) < bestErr) {
        bestErr = abs(fx);
        bestY = x;
      }
      if (prevF * fx <= 0.0) {
        a = prevX;
        b = x;
        fa = prevF;
        fb = fx;
        break;
      }
      prevX = x;
      prevF = fx;
    }
  }
  if (fa * fb > 0.0) { return bestY; }
  for (int it = 0; it < 60; ++it) {
    const double c = 0.5 * (a + b);
    const double fc = f(c);
    if (fa * fc <= 0.0) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }
  }
  return 0.5 * (a + b);
}

// -----------------------------------------------------------------
// AdrFixedF0 class
// -----------------------------------------------------------------

void QstlsF0Util::AdrFixedF0::get(const vector<double> &wvg, Vector3D &res) const {
  const int nx = wvg.size();
  const int nl = res.size(1);
  if (x == 0.0) { res.fill(0, 0.0); }
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
      auto func2 = [&](const double &t) -> double { return integrand2(t, wvg[i], l); };
      try {
        itg->compute(func1, func2, Itg2DParam(qMin, qMax, tMin, tMax), itgGrid);
        res(ix, l, i) = itg->getSolution();
      } catch (const std::exception &) {
        res(ix, l, i) = 0.0;
      }
    }
  }
}

double QstlsF0Util::AdrFixedF0::integrand1(const double &q, const double &l) const {
  if (l == 0) {
    if (q < 1.0e-8) { return 0.0; }
    return -0.5 * Theta * f0.derivative(q);
  }
  return q * f0.value(q);
}

double QstlsF0Util::AdrFixedF0::integrand2(const double &t,
                                           const double &y,
                                           const double &l) const {
  const double q = itg->getX();
  if (y == 0.0) { return 0.0; }
  const double x2 = x * x;
  const double y2 = y * y;
  const double q2 = q * q;
  const double txq = 2.0 * x * q;
  if (l == 0) {
    if (t == txq) { return 2.0 * q2 / (y2 + 2.0 * txq - x2); }
    if (x == y && t == 0.0) { return q; }
    const double t2 = t * t;
    double logarg = (t + txq) / (t - txq);
    logarg = (logarg < 0.0) ? -logarg : logarg;
    return y / (2.0 * t + y2 - x2)
           * ((q2 - t2 / (4.0 * x2)) * log(logarg) + q * t / x);
  }
  if (x == y && t == 0.0) { return 0.0; }
  const double tplT = 2.0 * M_PI * l * Theta;
  const double tplT2 = tplT * tplT;
  const double txqpt = txq + t;
  const double txqmt = txq - t;
  const double txqpt2 = txqpt * txqpt;
  const double txqmt2 = txqmt * txqmt;
  const double logarg = (txqpt2 + tplT2) / (txqmt2 + tplT2);
  return y / (2.0 * t + y * y - x * x) * log(logarg);
}
