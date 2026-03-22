#include "vs/vsbase.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;
using namespace GridPoints;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;

// -----------------------------------------------------------------
// VSBase
// -----------------------------------------------------------------

int VSBase::compute() {
  try {
    println("Free parameter calculation ...");
    doIterations();
    println("Done");
    return 0;
  } catch (const runtime_error &err) {
    cerr << err.what() << endl;
    return 1;
  }
}

const vector<vector<double>> &VSBase::getFreeEnergyIntegrand() const {
  return fxcIntegrand;
}

const vector<double> &VSBase::getFreeEnergyGrid() const { return rsGrid; }

void VSBase::doIterations() {
  auto func = [this](const double &alphaTmp) -> double {
    return alphaDifference(alphaTmp);
  };
  SecantSolver rsol(in().getErrMinAlpha(), in().getNIterAlpha());
  rsol.solve(func, in().getAlphaGuess());
  alpha = rsol.getSolution();
  println(formatUtil::format("Free parameter = {:.5f}", alpha));
}

double VSBase::alphaDifference(const double &alphaTmp) {
  alpha = alphaTmp;
  runGrid();
  updateFxcIntegrand();
  const double alphaTheoretical = computeAlpha();
  return alpha - alphaTheoretical;
}

double VSBase::computeAlpha() {
  const auto fed = getFreeEnergyData();
  const auto qdat = computeQData();
  const double fxcr = fed[1];
  const double fxcrr = fed[2];
  const double fxct = fed[3];
  const double fxctt = fed[4];
  const double fxcrt = fed[5];
  const double Q = qdat[0];
  const double Qr = qdat[1];
  const double Qt = qdat[2];
  const double numer = Q + (1.0 / 3.0) * fxcr - (1.0 / 6.0) * fxcrr
                       - (2.0 / 3.0) * (fxctt + fxcrt) + (1.0 / 3.0) * fxct;
  const double denom = Q + (1.0 / 3.0) * Qr + (2.0 / 3.0) * Qt;
  return numer / denom;
}

std::vector<double> VSBase::computeQData() {
  const double rs = getCoupling(GridPoints::CENTER);
  const double q = computeQRaw(GridPoints::CENTER) / rs;
  const double drs = getCoupling(GridPoints::RS_UP_THETA) - rs;
  const double qr = (computeQRaw(GridPoints::RS_UP_THETA)
                     - computeQRaw(GridPoints::RS_DOWN_THETA))
                        / (2.0 * drs)
                    - q;
  const double theta = getDegeneracy(GridPoints::CENTER);
  const double dt = getDegeneracy(GridPoints::RS_THETA_UP) - theta;
  const double qt = theta
                    * (computeQRaw(GridPoints::RS_THETA_UP) / rs
                       - computeQRaw(GridPoints::RS_THETA_DOWN) / rs)
                    / (2.0 * dt);
  return {q, qr, qt};
}

// -----------------------------------------------------------------
// rs grid and free energy integrand (from ThermoPropBase)
// -----------------------------------------------------------------

void VSBase::setRsGrid() {
  const double &rs = inScheme().getCoupling();
  const double &drs = in().getCouplingResolution();
  // rs must be a multiple of drs
  if (!numUtil::isZero(remainder(rs, drs))) {
    MPIUtil::throwError(
        "Inconsistent input parameters: the coupling parameter must be a "
        "multiple of the coupling resolution");
  }
  rsGrid.push_back(0.0);
  const double rsMax = rs + drs;
  while (!numUtil::equalTol(rsGrid.back(), rsMax)) {
    rsGrid.push_back(rsGrid.back() + drs);
  }
}

void VSBase::setFxcIntegrand() {
  const size_t nrs = rsGrid.size();
  fxcIntegrand.resize(3);
  for (auto &f : fxcIntegrand) {
    f.resize(nrs);
    vecUtil::fill(f, numUtil::Inf);
  }
  // Load pre-computed integrand from input if provided
  const auto &fxciData = in().getFreeEnergyIntegrand();
  if (!fxciData.grid.empty()) {
    for (int t = 0; t < 3; ++t) {
      const double rsMaxi = fxciData.grid.back();
      const Interpolator1D itp(fxciData.grid, fxciData.integrand[t]);
      if (itp.isValid()) {
        for (size_t i = 0; i < nrs; ++i) {
          if (rsGrid[i] <= rsMaxi) { fxcIntegrand[t][i] = itp.eval(rsGrid[i]); }
        }
      }
    }
  }
}

void VSBase::updateFxcIntegrand() {
  // Find the index of the target state point in rsGrid
  const double targetRs = getCoupling(CENTER);
  auto it = find_if(rsGrid.begin(), rsGrid.end(), [&](const double &rs) {
    return numUtil::equalTol(rs, targetRs);
  });
  if (it == rsGrid.end()) {
    MPIUtil::throwError(
        "Failed to find the target state point in the free energy grid");
  }
  const size_t idx = static_cast<size_t>(distance(rsGrid.begin(), it));
  // theta index 0 = DOWN, 1 = CENTER, 2 = UP
  fxcIntegrand[0][idx - 1] = getFxcIntegrandValue(RS_DOWN_THETA_DOWN);
  fxcIntegrand[0][idx] = getFxcIntegrandValue(RS_THETA_DOWN);
  fxcIntegrand[0][idx + 1] = getFxcIntegrandValue(RS_UP_THETA_DOWN);
  fxcIntegrand[1][idx - 1] = getFxcIntegrandValue(RS_DOWN_THETA);
  fxcIntegrand[1][idx] = getFxcIntegrandValue(CENTER);
  fxcIntegrand[1][idx + 1] = getFxcIntegrandValue(RS_UP_THETA);
  fxcIntegrand[2][idx - 1] = getFxcIntegrandValue(RS_DOWN_THETA_UP);
  fxcIntegrand[2][idx] = getFxcIntegrandValue(RS_THETA_UP);
  fxcIntegrand[2][idx + 1] = getFxcIntegrandValue(RS_UP_THETA_UP);
}

double VSBase::computeFreeEnergy(GridPoint p, bool normalize) const {
  // Map grid point to theta index (0=down, 1=center, 2=up)
  const int thetaIdx = static_cast<int>(p.theta) + 1;
  return thermoUtil::computeFreeEnergy(
      rsGrid, fxcIntegrand[thetaIdx], getCoupling(p), normalize);
}

vector<double> VSBase::getFreeEnergyData() const {
  // Free energy at the target state point
  const double fxc = computeFreeEnergy(CENTER, true);
  const double rs = getCoupling(CENTER);
  const double drs = getCoupling(RS_UP_THETA) - rs;
  // Derivatives w.r.t. coupling parameter
  double fxcr, fxcrr;
  {
    const double f0 = computeFreeEnergy(RS_UP_THETA, false);
    const double f1 = computeFreeEnergy(CENTER, false);
    const double f2 = computeFreeEnergy(RS_DOWN_THETA, false);
    fxcr = (f0 - f2) / (2.0 * drs * rs) - 2.0 * fxc;
    fxcrr = (f0 - 2.0 * f1 + f2) / (drs * drs) - 2.0 * fxc - 4.0 * fxcr;
  }
  // Derivatives w.r.t. degeneracy parameter
  double fxct, fxctt;
  {
    const double theta = getDegeneracy(CENTER);
    const double theta2 = theta * theta;
    const double dt = getDegeneracy(RS_THETA_UP) - theta;
    const double f0 = computeFreeEnergy(RS_THETA_UP, true);
    const double f1 = computeFreeEnergy(RS_THETA_DOWN, true);
    fxct = theta * (f0 - f1) / (2.0 * dt);
    fxctt = theta2 * (f0 - 2.0 * fxc + f1) / (dt * dt);
  }
  // Mixed derivative
  double fxcrt;
  {
    const double theta = getDegeneracy(CENTER);
    const double t_rs = theta / rs;
    const double dt = getDegeneracy(RS_THETA_UP) - theta;
    const double f0 = computeFreeEnergy(RS_UP_THETA_UP, false);
    const double f1 = computeFreeEnergy(RS_UP_THETA_DOWN, false);
    const double f2 = computeFreeEnergy(RS_DOWN_THETA_UP, false);
    const double f3 = computeFreeEnergy(RS_DOWN_THETA_DOWN, false);
    fxcrt = t_rs * (f0 - f1 - f2 + f3) / (4.0 * drs * dt) - 2.0 * fxct;
  }
  return {fxc, fxcr, fxcrr, fxct, fxctt, fxcrt};
}

GridPoint VSBase::getOutputGridPoint() const {
  const bool zeroCoupling = (inScheme().getCoupling() == 0.0);
  const bool zeroDegeneracy = (inScheme().getDegeneracy() == 0.0);
  if (zeroCoupling && zeroDegeneracy) return RS_DOWN_THETA_DOWN;
  if (!zeroCoupling && zeroDegeneracy) return RS_THETA_DOWN;
  if (zeroCoupling && !zeroDegeneracy) return RS_DOWN_THETA;
  return CENTER;
}

// -----------------------------------------------------------------
// QAdder
// -----------------------------------------------------------------

QAdder QAdder::classical(const std::vector<double> &wvg,
                         const std::vector<double> &ssf,
                         const std::shared_ptr<const Input> &in) {
  QAdder q;
  q.mode = Mode::CLASSICAL;
  q.classWvg = wvg;
  q.classSsf = ssf;
  q.classIn = in;
  return q;
}

QAdder QAdder::quantum(double Theta,
                       double mu,
                       double limitMin,
                       double limitMax,
                       const std::vector<double> &itgGrid,
                       std::shared_ptr<Integrator1D> itg1,
                       std::shared_ptr<Integrator2D> itg2,
                       std::shared_ptr<Interpolator1D> interp) {
  QAdder q;
  q.mode = Mode::QUANTUM;
  q.Theta = Theta;
  q.mu = mu;
  q.limits = {limitMin, limitMax};
  q.itgGridPtr = &itgGrid;
  q.itg1 = itg1;
  q.itg2 = itg2;
  q.interp = interp;
  return q;
}

double QAdder::get() const {
  if (mode == Mode::CLASSICAL) {
    return thermoUtil::computeInternalEnergy(
        classWvg, classSsf, 1.0, classIn->getDimension());
  }
  double Denominator;
  getIntDenominator(Denominator);
  auto func1 = [&](const double &w) -> double {
    return integrandNumerator2(w);
  };
  auto func2 = [&](const double &q) -> double {
    return integrandNumerator1(q);
  };
  itg2->compute(
      func1,
      func2,
      Itg2DParam(limits.first, limits.second, limits.first, limits.second),
      *itgGridPtr);
  return 12.0 / (M_PI * numUtil::lambda) * itg2->getSolution() / Denominator;
}

double QAdder::ssf(const double &y) const { return interp->eval(y); }

double QAdder::integrandDenominator(const double y) const {
  const double y2 = y * y;
  return 1.0 / (exp(y2 / Theta - mu) + 1.0);
}

double QAdder::integrandNumerator1(const double q) const {
  const double w = itg2->getX();
  if (q == 0.0) { return 0.0; }
  const double w2 = w * w;
  const double q2 = q * q;
  double logarg = (w + 2 * q) / (w - 2 * q);
  logarg = (logarg < 0.0) ? -logarg : logarg;
  if (w == 0.0) { return 1.0 / (12.0 * (exp(q2 / Theta - mu) + 1.0)); }
  return q2 / (exp(q2 / Theta - mu) + 1.0) * (q / w * log(logarg) - 1.0) / w2;
}

double QAdder::integrandNumerator2(const double w) const {
  return (ssf(w) - 1.0);
}

void QAdder::getIntDenominator(double &res) const {
  auto func = [&](double y) -> double { return integrandDenominator(y); };
  itg1->compute(func, ItgParam(limits.first, limits.second));
  res = itg1->getSolution();
}
