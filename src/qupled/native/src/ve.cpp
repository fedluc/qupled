#include "ve.hpp"
#include "format.hpp"
#include "mpi_util.hpp"
#include "thermo_util.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

using namespace std;
using namespace MPIUtil;
using Itg2DParam = Integrator2D::Param;

namespace {

constexpr double EPS = 1.0e-12;

bool samePoint(const double a, const double b) {
  return std::abs(a - b) <= 1.0e-10;
}

double average(const std::vector<double> &values) {
  if (values.empty()) { return 0.0; }
  double sum = 0.0;
  for (const double value : values) { sum += value; }
  return sum / static_cast<double>(values.size());
}

} // namespace

// -----------------------------------------------------------------
// VE class
// -----------------------------------------------------------------

VE::VE(const std::shared_ptr<const VEInput> &in_)
    : Qstls(in_, true),
      a0Coeff(0.0),
      kxc0(0.0),
      mxcInf(0.0),
      cxc(0.0),
      omegaM2(1.0) {
  if (in().getDimension() == dimensionsUtil::Dimension::D2) {
    throwError("2D calculations are not implemented for this scheme.");
  }
  if (in().getDegeneracy() <= 0.0) {
    throwError("VE is implemented only for finite degeneracy (theta > 0)");
  }
  const size_t nx = wvg.size();
  const size_t nl = in().getNMatsubara();
  adrFixedVe.resize(nx, nl, nx);
  aCoeff.assign(nl, 0.0);
  a1Coeff.assign(nl, 0.0);
  mxcLDiag.assign(nl, 0.0);
  matsubaraGrid.resize(nl);
  for (size_t l = 0; l < nl; ++l) { matsubaraGrid[l] = static_cast<double>(l); }
}

void VE::init() {
  Qstls::init();
  print("Computing viscoelastic fixed kernel: ");
  fflush(stdout);
  computeAdrFixedVe();
  println("Done");
  initializeFreeEnergyIntegrand();
}

void VE::initializeFreeEnergyIntegrand() {
  const auto &fxcInput = in().getFreeEnergyIntegrand();
  freeEnergyGridOut = fxcInput.grid;
  freeEnergyIntegrandOut = fxcInput.integrand;
  if (freeEnergyIntegrandOut.empty()) { freeEnergyIntegrandOut.resize(1); }
  const double hfZero = exactZeroCouplingIntegrand();
  if (freeEnergyGridOut.empty() || !samePoint(freeEnergyGridOut.front(), 0.0)) {
    freeEnergyGridOut.insert(freeEnergyGridOut.begin(), 0.0);
    freeEnergyIntegrandOut[0].insert(freeEnergyIntegrandOut[0].begin(), hfZero);
  } else if (!freeEnergyIntegrandOut[0].empty()) {
    freeEnergyIntegrandOut[0][0] = hfZero;
  }
  if (freeEnergyIntegrandOut[0].size() != freeEnergyGridOut.size()) {
    freeEnergyIntegrandOut[0].assign(freeEnergyGridOut.size(), 0.0);
    freeEnergyIntegrandOut[0][0] = hfZero;
  }
  updateFreeEnergyIntegrandCurrent();
}

void VE::updateFreeEnergyIntegrandCurrent() {
  if (freeEnergyIntegrandOut.empty()) { freeEnergyIntegrandOut.resize(1); }
  if (freeEnergyIntegrandOut[0].size() != freeEnergyGridOut.size()) {
    freeEnergyIntegrandOut[0].resize(freeEnergyGridOut.size(), 0.0);
  }
  const double rs = in().getCoupling();
  const double value = currentFreeEnergyIntegrand();
  auto it = std::lower_bound(freeEnergyGridOut.begin(), freeEnergyGridOut.end(), rs);
  if (it != freeEnergyGridOut.end() && samePoint(*it, rs)) {
    const size_t idx = distance(freeEnergyGridOut.begin(), it);
    freeEnergyIntegrandOut[0][idx] = value;
    return;
  }
  const size_t idx = distance(freeEnergyGridOut.begin(), it);
  freeEnergyGridOut.insert(it, rs);
  freeEnergyIntegrandOut[0].insert(freeEnergyIntegrandOut[0].begin() + idx, value);
}

double VE::currentFreeEnergyIntegrand() const {
  return thermoUtil::computeInternalEnergy(wvg, ssf, 1.0, in().getDimension());
}

double VE::exactZeroCouplingIntegrand() const {
  return thermoUtil::computeInternalEnergy(wvg, ssfHF, 1.0, in().getDimension());
}

double VE::computeStaticKxc0() const {
  const auto &grid = freeEnergyGridOut;
  const auto &integrand = freeEnergyIntegrandOut[0];
  const size_t n = grid.size();
  if (n < 3) { return 0.0; }

  const double rs = in().getCoupling();
  auto it = std::lower_bound(grid.begin(), grid.end(), rs);
  if (it == grid.end() || !samePoint(*it, rs)) { return 0.0; }
  const size_t idx = distance(grid.begin(), it);

  auto freeEnergy = [&](size_t i) {
    const double coupling = grid[i];
    if (samePoint(coupling, 0.0)) {
      return 0.5 * integrand[i];
    }
    return thermoUtil::computeFreeEnergy(grid, integrand, coupling, true);
  };

  double fp = 0.0;
  double fpp = 0.0;

  if (idx >= 1 && idx + 1 < n) {
    const double h = grid[idx + 1] - grid[idx];
    const double f0 = freeEnergy(idx + 1);
    const double f1 = freeEnergy(idx);
    const double f2 = freeEnergy(idx - 1);
    fp = (f0 - f2) / (2.0 * h);
    fpp = (f0 - 2.0 * f1 + f2) / (h * h);
  } else if (idx >= 2) {
    const double h = grid[idx] - grid[idx - 1];
    const double f0 = freeEnergy(idx);
    const double f1 = freeEnergy(idx - 1);
    const double f2 = freeEnergy(idx - 2);
    fp = (3.0 * f0 - 4.0 * f1 + f2) / (2.0 * h);
    fpp = (f0 - 2.0 * f1 + f2) / (h * h);
  } else if (idx + 2 < n) {
    const double h = grid[idx + 1] - grid[idx];
    const double f0 = freeEnergy(idx);
    const double f1 = freeEnergy(idx + 1);
    const double f2 = freeEnergy(idx + 2);
    fp = (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h);
    fpp = (f0 - 2.0 * f1 + f2) / (h * h);
  } else {
    return 0.0;
  }

  return (rs * rs * fpp - 2.0 * rs * fp) / 9.0;
}

double VE::computeSmallXCoefficient(const Vector2D &psi, const size_t l) const {
  double numer = 0.0;
  double denom = 0.0;
  const size_t nFit = std::min<size_t>(5, wvg.size() > 0 ? wvg.size() - 1 : 0);
  for (size_t i = 1; i <= nFit; ++i) {
    const double x = wvg[i];
    const double phi = idr(i, l);
    if (std::abs(phi) < EPS) { continue; }
    const double x2 = x * x;
    const double ratio = psi(i, l) / phi;
    numer += x2 * ratio;
    denom += x2 * x2;
  }
  return (std::abs(denom) < EPS) ? 0.0 : numer / denom;
}

std::vector<double> VE::computeA1(const Vector2D &psi1) const {
  const size_t nl = psi1.size(1);
  std::vector<double> a1(nl, 0.0);
  for (size_t l = 0; l < nl; ++l) { a1[l] = computeSmallXCoefficient(psi1, l); }
  return a1;
}

double VE::estimateCxcInitial(const std::vector<double> &a1) const {
  const double rs = in().getCoupling();
  const double theta = in().getDegeneracy();
  const size_t nl = a1.size();
  if (nl <= 1) { return 0.0; }
  const size_t start = std::max<size_t>(1, (2 * nl) / 3);
  std::vector<double> plateau;
  plateau.reserve(nl - start);
  for (size_t l = start; l < nl; ++l) {
    const double u = 2.0 * M_PI * static_cast<double>(l) * theta;
    const double value = -(3.0 * numUtil::lambda * numUtil::lambda / rs) * u * u * a1[l];
    if (std::isfinite(value)) { plateau.push_back(std::abs(value)); }
  }
  return average(plateau);
}

double VE::estimateCxcFromMxc(const std::vector<double> &mxcL) const {
  const double theta = in().getDegeneracy();
  const size_t nl = mxcL.size();
  if (nl <= 1) { return 0.0; }
  const size_t start = std::max<size_t>(1, (2 * nl) / 3);
  std::vector<double> plateau;
  plateau.reserve(nl - start);
  for (size_t l = start; l < nl; ++l) {
    const double u = 2.0 * M_PI * static_cast<double>(l) * theta;
    const double value = u * u * (mxcInf - mxcL[l]);
    if (std::isfinite(value)) { plateau.push_back(value); }
  }
  return average(plateau);
}

void VE::computeLfc() {
  if (in().getDegeneracy() == 0.0) { return; }
  const int nx = static_cast<int>(wvg.size());
  const int nl = static_cast<int>(in().getNMatsubara());
  const auto ssfi = make_shared<Interpolator1D>(wvg, ssf);
  Vector2D psi0(nx, nl);
  Vector2D psi1(nx, nl);

  for (int i = 0; i < nx; ++i) {
    QstlsUtil::Adr adr0(in().getDegeneracy(), wvg.front(), wvg.back(), wvg[i], ssfi, itg);
    adr0.get(wvg, adrFixed, psi0);
    QstlsUtil::Adr adr1(in().getDegeneracy(), wvg.front(), wvg.back(), wvg[i], ssfi, itg);
    adr1.get(wvg, adrFixedVe, psi1);
  }

  updateFreeEnergyIntegrandCurrent();

  const double rs = in().getCoupling();
  const double theta = in().getDegeneracy();
  const double uInt = getUInt();
  const double a0 = -0.5 * M_PI * numUtil::lambda * rs * uInt;
  a0Coeff = a0;
  kxc0 = computeStaticKxc0();
  mxcInf = 1.5 * M_PI * std::pow(numUtil::lambda, 3) * uInt;

  const std::vector<double> a1 = computeA1(psi1);
  a1Coeff = a1;
  const double deltaM = mxcInf - kxc0;
  cxc = estimateCxcInitial(a1);
  if (!std::isfinite(cxc) || std::abs(cxc) < EPS) { cxc = 1.0; }

  std::vector<double> mxcL(nl, kxc0);
  constexpr int maxInnerIter = 30;
  constexpr double cxcTol = 1.0e-4;
  for (int iter = 0; iter < maxInnerIter; ++iter) {
    omegaM2 = (std::abs(deltaM) < EPS) ? 1.0 : std::abs(cxc / deltaM);
    for (int l = 0; l < nl; ++l) {
      const double u = 2.0 * M_PI * static_cast<double>(l) * theta;
      const double rl = (l == 0) ? 0.0 : (u * u) / (u * u + omegaM2);
      const double mTarget = kxc0 + deltaM * rl;
      mxcL[l] = mTarget;
      const double aTarget =
          -(rs / (3.0 * numUtil::lambda * numUtil::lambda)) * mTarget;
      if (std::abs(a1[l]) < EPS) {
        aCoeff[l] = 0.0;
      } else {
        aCoeff[l] = (aTarget - a0) / a1[l];
      }
    }

    for (int l = 0; l < nl; ++l) {
      const double aFull = a0 + aCoeff[l] * a1[l];
      mxcL[l] = -(3.0 * numUtil::lambda * numUtil::lambda / rs) * aFull;
    }

    const double updatedCxc = estimateCxcFromMxc(mxcL);
    if (!std::isfinite(updatedCxc)) { break; }
    const double denom = std::max(std::abs(cxc), 1.0);
    const double relErr = std::abs(updatedCxc - cxc) / denom;
    cxc = 0.5 * cxc + 0.5 * updatedCxc;
    if (relErr < cxcTol) { break; }
  }
  mxcLDiag = mxcL;

  lfc = psi0;
  for (int i = 0; i < nx; ++i) {
    for (int l = 0; l < nl; ++l) {
      lfc(i, l) += aCoeff[l] * psi1(i, l);
    }
  }
  lfc.div(idr);
  lfc.fill(0, 0.0);
}

void VE::computeAdrFixedVe() {
  if (in().getDegeneracy() == 0.0) { return; }
  const int nx = wvg.size();
  const int nl = in().getNMatsubara();
  const string name = formatUtil::format("{}_VE1", adrFixedDatabaseName);
  if (in().getFixedRunId() != DEFAULT_INT) {
    adrFixedVe.resize(nx, nl, nx);
    readAdrFixed(adrFixedVe, name, in().getFixedRunId());
    return;
  }

  const int nxnl = nx * nl;
  const bool segregatedItg = in().getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  atomic<int> completed{0};
  const int progressStep = max(1, nx / 10);

  auto loopFunc = [&](int i) -> void {
    shared_ptr<Integrator2D> itg2 = make_shared<Integrator2D>(in().getIntError());
    VEUtil::AdrFixed adrTmp(in().getDegeneracy(),
                            wvg.front(),
                            wvg.back(),
                            wvg[i],
                            mu,
                            itgGrid,
                            itg2);
    adrTmp.get(wvg, adrFixedVe);
    const int done = ++completed;
    if (isRoot() && (done % progressStep == 0 || done == nx)) {
      cerr << "[VE] fixed kernel progress: " << done << "/" << nx << endl;
    }
  };

  const auto &loopData = parallelFor(loopFunc, nx, in().getNThreads());
  gatherLoopData(adrFixedVe.data(), loopData, nxnl);
  if (isRoot()) { writeAdrFixed(adrFixedVe, name); }
}

// -----------------------------------------------------------------
// VE fixed kernel
// -----------------------------------------------------------------

void VEUtil::AdrFixed::get(const vector<double> &wvg, Vector3D &res) const {
  const int nx = wvg.size();
  const int nl = res.size(1);
  if (x == 0.0) {
    auto it0 = find(wvg.begin(), wvg.end(), x);
    if (it0 != wvg.end()) { res.fill(distance(wvg.begin(), it0), 0.0); }
    return;
  }
  auto it = find(wvg.begin(), wvg.end(), x);
  assert(it != wvg.end());
  const size_t ix = distance(wvg.begin(), it);
  for (int l = 0; l < nl; ++l) {
    for (int i = 0; i < nx; ++i) {
      const double w = wvg[i];
      if (w == 0.0) {
        res(ix, l, i) = 0.0;
        continue;
      }
      auto func1 = [&](const double &y) -> double { return integrand1(y, l); };
      auto func2 = [&](const double &s) -> double { return integrand2(s, w, l); };
      try {
        itg->compute(func1, func2, Itg2DParam(qMin, qMax, -1.0, 1.0), itgGrid);
        res(ix, l, i) = itg->getSolution();
      } catch (const std::exception &) {
        res(ix, l, i) = 0.0;
      }
    }
  }
}

double VEUtil::AdrFixed::integrand1(const double &y, const double &l) const {
  const double expArg = std::exp(y * y / Theta - mu);
  if (l == 0) {
    const double denom = expArg + 1.0;
    return y * expArg / (denom * denom);
  }
  return y / (expArg + 1.0);
}

double VEUtil::AdrFixed::q(const double &s) const {
  return 5.0 / 9.0 + 4.0 * s * s / 3.0;
}

double VEUtil::AdrFixed::integrand2(const double &s,
                                    const double &w,
                                    const double &l) const {
  const double y = itg->getX();
  if (x == 0.0 || w == 0.0) { return 0.0; }
  const double x2 = x * x;
  const double xw = x * w;
  const double xq = x2 + xw * s;
  const double denom = w * w + x2 + 2.0 * xw * s;
  if (std::abs(denom) < EPS) { return 0.0; }
  const double pref = xw / denom;
  if (l == 0) {
    const double num1 = xq + 2.0 * x * y;
    const double num2 = xq - 2.0 * x * y;
    if (std::abs(num2) < EPS) { return 0.0; }
    double logarg = num1 / num2;
    logarg = (logarg < 0.0) ? -logarg : logarg;
    if (logarg <= 0.0) { return 0.0; }
    const double term1 = (y * y - xq * xq / (4.0 * x2)) * std::log(logarg);
    const double term2 = y * xq / x;
    return pref * q(s) * (term1 + term2);
  }
  const double u = 2.0 * M_PI * static_cast<double>(l) * Theta;
  const double arg1 = 2.0 * x * y + xq;
  const double arg2 = 2.0 * x * y - xq;
  const double logarg = (arg1 * arg1 + u * u) / (arg2 * arg2 + u * u);
  return pref * q(s) * std::log(logarg);
}
