#include "qstlspimc.hpp"
#include "esa.hpp"
#include "mpi_util.hpp"
#include "num_util.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

using namespace std;
using namespace MPIUtil;
using Itg2DParam = Integrator2D::Param;

namespace {

  constexpr double TOL = 1.0e-8;
  constexpr double DEFAULT_A_CUTOFF = 0.5;

  bool eqTol(const double a, const double b, const double tol = TOL) {
    return abs(a - b) <= tol;
  }

  bool isSupportedRs(const double rs) {
    return eqTol(rs, 3.23, 1.0e-6) || eqTol(rs, 10.0) || eqTol(rs, 20.0);
  }

  string selectMomentumDistributionFile(const double rs) {
    if (eqTol(rs, 3.23, 1.0e-6)) { return "Momentum_UEG_rs3p23_theta1.txt"; }
    if (eqTol(rs, 10.0)) { return "Momentum_UEG_rs10_theta1.txt"; }
    if (eqTol(rs, 20.0)) { return "Momentum_UEG_rs20_theta1.txt"; }
    throwError("QSTLS-PIMC supports only rs = 3.23, 10 or 20 at theta = 1");
    return "";
  }

  vector<string> candidateDataPaths(const string &filename) {
    vector<string> paths;
    if (const char *env = std::getenv("QUPLED_PIMC_DATA_DIR")) {
      paths.push_back((filesystem::path(env) / filename).string());
    }
    paths.push_back(filename);
    paths.push_back((filesystem::path("src") / filename).string());
    paths.push_back((filesystem::path("..") / filename).string());
    const filesystem::path srcDir = filesystem::path(__FILE__).parent_path();
    paths.push_back((srcDir / filename).string());
    paths.push_back((srcDir / ".." / ".." / ".." / ".." / filename).string());
    paths.push_back((srcDir / ".." / ".." / ".." / ".." / "src" / filename)
                        .string());
    return paths;
  }

  pair<vector<double>, vector<double>> loadMomentumDistribution(const string &path) {
    ifstream in(path);
    if (!in) { throwError("Failed to open PIMC momentum distribution file: " + path); }
    vector<double> y;
    vector<double> f0;
    vector<double> yAll;
    vector<double> f0All;
    string line;
    while (getline(in, line)) {
      if (line.empty() || line[0] == '#') { continue; }
      istringstream iss(line);
      double yi = 0.0;
      double fi = 0.0;
      double yLow = 0.0;
      double yHigh = 0.0;
      string kind;
      if (!(iss >> yi >> fi)) { continue; }
      // Keep non-negative grid values and clip tiny negative noise in f0.
      if (yi < 0.0) { continue; }
      yAll.push_back(yi);
      f0All.push_back(max(fi, 0.0));
      if ((iss >> yLow >> yHigh >> kind)) {
        // Use only reliable points inside the Fermi surface when tagged.
        if (kind == "i" || kind == "I") {
          y.push_back(yi);
          f0.push_back(max(fi, 0.0));
        }
      } else {
        y.push_back(yi);
        f0.push_back(max(fi, 0.0));
      }
    }
    if (y.size() < 4) {
      y = yAll;
      f0 = f0All;
    }
    if (y.size() < 4) {
      throwError("Insufficient PIMC momentum distribution points in " + path);
    }
    return {y, f0};
  }

  string findMomentumDistributionFile(const string &filename) {
    for (const auto &path : candidateDataPaths(filename)) {
      if (filesystem::exists(path)) { return path; }
    }
    throwError("Unable to locate PIMC momentum distribution file " + filename
               + ". Expected in project root or src/.");
    return "";
  }

} // namespace

// -----------------------------------------------------------------
// QSTLS-PIMC class
// -----------------------------------------------------------------

QstlsPimc::QstlsPimc(const std::shared_ptr<const QstlsPimcInput> &in_)
    : Qstls(in_, true) {
  const double theta = in().getDegeneracy();
  const double rs = in().getCoupling();
  if (!eqTol(theta, 1.0) || !isSupportedRs(rs)) {
    throwError(
        "QSTLS-PIMC is implemented only for theta = 1 and rs = 3.23, 10 or 20");
  }
}

void QstlsPimc::init() {
  Stls::init();
  print("Computing fixed component of the auxiliary density response (PIMC): ");
  fflush(stdout);
  computeAdrFixedPimc();
  println("Done");
}

void QstlsPimc::computeAdrFixedPimc() {
  if (in().getDegeneracy() == 0.0) { return; }
  const int nx = wvg.size();
  const int nl = in().getNMatsubara();
  const QstlsPimcUtil::Distribution f0(in().getCoupling(),
                                       in().getDegeneracy(),
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
    QstlsPimcUtil::AdrFixedPimc adrTmp(in().getDegeneracy(),
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
      cerr << "[QSTLS-PIMC] fixed ADR progress: " << done << "/" << nx << endl;
    }
  };
  const auto &loopData = parallelFor(loopFunc, nx, in().getNThreads());
  gatherLoopData(adrFixed.data(), loopData, nxnl);
  if (isRoot()) { writeAdrFixed(adrFixed, adrFixedDatabaseName); }
}

void QstlsPimc::computeSsfFinite() {
  // PIMC-based qSTLS update:
  // S(x) = (3/2) * Theta * sum_l w_l * Phi(x,l) / [1 + v(x) * (Phi(x,l)-Psi(x,l))]
  // with Phi -> idr and Psi/Phi -> lfc in the current implementation.
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
// Momentum distribution helper
// -----------------------------------------------------------------

QstlsPimcUtil::Distribution::Distribution(const double &rs_,
                                          const double &theta_,
                                          const double &etaIn,
                                          const double &ySecIn,
                                          const double &aCutoffIn)
    : rs(rs_),
      theta(theta_),
      yGrid(),
      fGrid(),
      dfGrid(),
      yMin(0.0),
      yMax(0.0),
      yJoin(0.0),
      eta(0.0),
      ySec(0.0),
      aCutoff(0.0),
      tailCoeff(0.0) {
  if (!eqTol(theta, 1.0) || !isSupportedRs(rs)) {
    throwError("QSTLS-PIMC supports only theta = 1 and rs = 3.23, 10 or 20");
  }
  const string filename = selectMomentumDistributionFile(rs);
  const string path = findMomentumDistributionFile(filename);
  const auto [y, f0] = loadMomentumDistribution(path);
  yGrid = y;
  fGrid = f0;
  yMin = y.front();
  yMax = y.back();
  yJoin = 0.0;
  // Use ESA parametrization for g(0, rs, theta) in the asymptotic tail.
  const double g0 = ESAUtil::Slfc(rs, theta).onTop();
  const double rs2 = rs * rs;
  const double lambda2 = numUtil::lambda * numUtil::lambda;
  tailCoeff = 8.0 / (9.0 * M_PI * M_PI) * lambda2 * rs2 * abs(g0);

  // Convert PIMC momentum axis to the y-normalization used for the tail:
  // y = lambda * rs * x
  const double axisScale = numUtil::lambda * rs;
  for (auto &x : yGrid) { x *= axisScale; }
  yMin = yGrid.front();
  yMax = yGrid.back();

  // Precompute derivative on the normalized grid using non-uniform finite differences.
  dfGrid.assign(fGrid.size(), 0.0);
  if (fGrid.size() >= 3) {
    for (size_t i = 1; i + 1 < fGrid.size(); ++i) {
      const double dxm = yGrid[i] - yGrid[i - 1];
      const double dxp = yGrid[i + 1] - yGrid[i];
      const double sm = (fGrid[i] - fGrid[i - 1]) / dxm;
      const double sp = (fGrid[i + 1] - fGrid[i]) / dxp;
      dfGrid[i] = (dxp * sm + dxm * sp) / (dxm + dxp);
    }
    dfGrid.front() = 0.0;
    dfGrid.back() =
        (fGrid.back() - fGrid[fGrid.size() - 2]) / (yGrid.back() - yGrid[yGrid.size() - 2]);
  }

  // Smooth crossover from tabulated PIMC data to asymptotic 1/y^8 behavior.
  const bool hasEta = !std::isnan(etaIn);
  const bool hasYSec = !std::isnan(ySecIn);
  const bool hasACutoff = !std::isnan(aCutoffIn);
  eta = hasEta ? etaIn : 8.0;
  aCutoff = hasACutoff ? aCutoffIn : DEFAULT_A_CUTOFF;
  ySec = hasYSec ? ySecIn : determineYSec();
  const double m0 = moment(eta, ySec);
  if (isRoot()) {
    cerr << "[QSTLS-PIMC] normalization moment: " << m0
         << " target=" << (1.0 / 3.0) << endl;
  }
}

double QstlsPimcUtil::Distribution::pimcValue(const double &y) const {
  const double yc = std::clamp(y, yMin, yMax);
  const auto it = lower_bound(yGrid.begin(), yGrid.end(), yc);
  if (it == yGrid.begin()) { return fGrid.front(); }
  if (it == yGrid.end()) { return fGrid.back(); }
  const size_t i1 = distance(yGrid.begin(), it);
  const size_t i0 = i1 - 1;
  const double x0 = yGrid[i0];
  const double x1 = yGrid[i1];
  const double t = (yc - x0) / (x1 - x0);
  return (1.0 - t) * fGrid[i0] + t * fGrid[i1];
}

double QstlsPimcUtil::Distribution::pimcDerivative(const double &y) const {
  const double yc = std::clamp(abs(y), yMin, yMax);
  if (yc <= yMin + 1.0e-12) { return 0.0; }
  const auto it = lower_bound(yGrid.begin(), yGrid.end(), yc);
  if (it == yGrid.begin()) { return 0.0; }
  if (it == yGrid.end()) { return dfGrid.back(); }
  const size_t i1 = distance(yGrid.begin(), it);
  const size_t i0 = i1 - 1;
  const double x0 = yGrid[i0];
  const double x1 = yGrid[i1];
  const double t = (yc - x0) / (x1 - x0);
  return (1.0 - t) * dfGrid[i0] + t * dfGrid[i1];
}

double QstlsPimcUtil::Distribution::asymptoticValue(const double &y) const {
  const double yc = max(y, 1.0e-8);
  const double y2 = yc * yc;
  const double y4 = y2 * y2;
  return tailCoeff / (y4 * y4);
}

double QstlsPimcUtil::Distribution::asymptoticDerivative(const double &y) const {
  const double yc = max(y, 1.0e-8);
  return -8.0 * tailCoeff / pow(yc, 9);
}

double QstlsPimcUtil::Distribution::switchFunction(const double &y) const {
  if (y <= aCutoff) { return 0.0; }
  const double sw = 0.5 * (1.0 + tanh(eta * (y - ySec)));
  return sw;
}

double QstlsPimcUtil::Distribution::switchDerivative(const double &y) const {
  if (y <= aCutoff) { return 0.0; }
  const double arg = eta * (y - ySec);
  const double c = cosh(arg);
  return 0.5 * eta / (c * c);
}

double QstlsPimcUtil::Distribution::value(const double &y) const {
  return value(y, eta, ySec);
}

double QstlsPimcUtil::Distribution::value(const double &y,
                                          const double &etaLocal,
                                          const double &ySecLocal) const {
  const double fp = pimcValue(y);
  const double fa = asymptoticValue(y);
  if (y <= aCutoff) { return fp; }
  const double a = 0.5 * (1.0 + tanh(etaLocal * (y - ySecLocal)));
  return (1.0 - a) * fp + a * fa;
}

double QstlsPimcUtil::Distribution::determineYSec() const {
  // Two-y strategy: yJoin is fixed, and ySec is solved from normalization.
  return solveYSec(eta);
}

double QstlsPimcUtil::Distribution::derivative(const double &y) const {
  const double fp = pimcValue(y);
  const double fa = asymptoticValue(y);
  const double dfp = pimcDerivative(y);
  const double dfa = asymptoticDerivative(y);
  const double a = switchFunction(y);
  const double da = switchDerivative(y);
  return (1.0 - a) * dfp + a * dfa + da * (fa - fp);
}

double QstlsPimcUtil::Distribution::moment(const double &etaLocal,
                                           const double &ySecLocal) const {
  // Normalization in reduced units: \int_0^\infty y^2 f0(y) dy = 1/3
  const double yUpper = max(100.0, 8.0 * yMax);
  const int n = 8000;
  const double h = yUpper / static_cast<double>(n);
  double sum = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double y = i * h;
    const double w = (i == 0 || i == n) ? 0.5 : 1.0;
    sum += w * y * y * value(y, etaLocal, ySecLocal);
  }
  // Add analytic tail above yUpper
  const double tail = tailCoeff / (5.0 * pow(yUpper, 5));
  return h * sum + tail;
}

double QstlsPimcUtil::Distribution::solveYSec(const double &etaLocal) const {
  const double target = 1.0 / 3.0;
  auto f = [&](double ySecLocal) { return moment(etaLocal, ySecLocal) - target; };
  double a = 0.01 * yMax;
  double b = 10.0 * yMax;
  double fa = f(a);
  double fb = f(b);
  double bestY = a;
  double bestErr = abs(fa);
  if (fa * fb > 0.0) {
    const int nScan = 220;
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
  if (fa * fb > 0.0) {
    return bestY;
  }
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

double QstlsPimcUtil::Distribution::solveEta(const double &ySecLocal) const {
  const double target = 1.0 / 3.0;
  auto f = [&](double etaLocal) { return moment(etaLocal, ySecLocal) - target; };
  double a = 0.1;
  double b = 200.0;
  double fa = f(a);
  double fb = f(b);
  if (fa * fb > 0.0) {
    throwError("Failed to bracket eta for PIMC normalization");
  }
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
// AdrFixedPimc class
// -----------------------------------------------------------------

void QstlsPimcUtil::AdrFixedPimc::get(const vector<double> &wvg, Vector3D &res) const {
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
      auto func2 = [&](const double &t) -> double {
        return integrand2(t, wvg[i], l);
      };
      try {
        itg->compute(func1, func2, Itg2DParam(qMin, qMax, tMin, tMax), itgGrid);
        res(ix, l, i) = itg->getSolution();
      } catch (const std::exception &) {
        // Regularize isolated divergent quadrature points.
        res(ix, l, i) = 0.0;
      }
    }
  }
}

double QstlsPimcUtil::AdrFixedPimc::integrand1(const double &q, const double &l) const {
  if (l == 0) {
    if (q < 1.0e-8) { return 0.0; }
    return -0.5 * Theta * f0.derivative(q);
  }
  return q * f0.value(q);
}

double QstlsPimcUtil::AdrFixedPimc::integrand2(const double &t,
                                               const double &y,
                                               const double &l) const {
  const double q = itg->getX();
  if (y == 0) { return 0.0; }
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
