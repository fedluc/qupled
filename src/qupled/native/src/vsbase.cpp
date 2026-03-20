#include "vsbase.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;
using namespace GridPoints;

// -----------------------------------------------------------------
// VSBase
// -----------------------------------------------------------------

int VSBase::compute() {
  try {
    init();
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
  updateSolution();
}

double VSBase::alphaDifference(const double &alphaTmp) {
  alpha = alphaTmp;
  runGrid();
  updateFxcIntegrand();
  const double alphaTheoretical = computeAlpha();
  return alpha - alphaTheoretical;
}

// Unified alpha formula. Works for both VSStls and QVSStls because the
// classical Q-term (internal energy) is a limiting case of the quantum one.
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
  // Free energy at the target state point (normalised)
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
