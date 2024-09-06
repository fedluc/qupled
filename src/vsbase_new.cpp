#include "vsbase_new.hpp"
#include "input.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;

// -----------------------------------------------------------------
// VSBase class
// -----------------------------------------------------------------

int VSBase::compute() {
  try {
    init();
    if (verbose) cout << "Free parameter calculation ..." << endl;
    doIterations();
    if (verbose) cout << "Done" << endl;
    return 0;
  } catch (const runtime_error &err) {
    cerr << err.what() << endl;
    return 1;
  }
}

void VSBase::init() {
  initScheme();
  initFreeEnergyIntegrand();
}

vector<vector<double>> VSBase::getFreeEnergyIntegrand() const {
  return getThermoProp().getFreeEnergyIntegrand();
}

vector<double> VSBase::getFreeEnergyGrid() const {
  return getThermoProp().getFreeEnergyGrid();
}

vector<double> VSBase::getAlpha() const {
  return getThermoProp().getAlpha();
}

void VSBase::doIterations() {
  auto func = [this](const double &alphaTmp) -> double {
    return alphaDifference(alphaTmp);
  };
  SecantSolver rsol(in.getErrMinAlpha(), in.getNIterAlpha());
  rsol.solve(func, in.getAlphaGuess());
  alpha = rsol.getSolution();
  if (verbose) { std::cout << "Free parameter = " << alpha << std::endl; }
  updateSolution();
}

double VSBase::alphaDifference(const double &alphaTmp) {
  alpha = alphaTmp;
  getThermoProp().setAlpha(alpha);
  const double alphaTheoretical = computeAlpha();
  return alpha - alphaTheoretical;
}

// -----------------------------------------------------------------
// ThermoPropBase class
// -----------------------------------------------------------------

ThermoPropBase::ThermoPropBase(const VSInput &in)
  : verbose(MPIUtil::isRoot()) { //, structProp(in) {
  // Check if we are solving for particular state points
  isZeroCoupling = (in.getCoupling() == 0.0);
  isZeroDegeneracy = (in.getDegeneracy() == 0.0);
  // Set variables related to the free energy calculation
  setRsGrid(in);
  setFxcIntegrand(in);
  setAlpha(in);
  setFxcIdxTargetStatePoint(in);
  setFxcIdxUnsolvedStatePoint();
}


void ThermoPropBase::setRsGrid(const VSInput &in) {
  const double &rs = in.getCoupling();
  const double &drs = in.getCouplingResolution();
  if (!numUtil::isZero(std::remainder(rs, drs))) {
    MPIUtil::throwError("Inconsistent input parameters: the coupling parameter must be a "
			"multiple of the coupling resolution");
  }
  rsGrid.push_back(0.0);
  const double rsMax = rs + drs;
  while (!numUtil::equalTol(rsGrid.back(), rsMax)) {
    rsGrid.push_back(rsGrid.back() + drs);
  }
}

void ThermoPropBase::setFxcIntegrand(const VSInput &in) {
  const size_t nrs = rsGrid.size();
  // Initialize the free energy integrand
  fxcIntegrand.resize(NPOINTS);
  for (auto &f : fxcIntegrand) {
    f.resize(nrs);
    vecUtil::fill(f, numUtil::Inf);
  }
  // Fill the free energy integrand and the free parameter if passed in input
  const auto &fxciData = in.getFreeEnergyIntegrand();
  if (!fxciData.grid.empty()) {
    for (const auto &theta : {Idx::THETA_DOWN, Idx::THETA, Idx::THETA_UP}) {
      const double rsMaxi = fxciData.grid.back();
      const Interpolator1D itp(fxciData.grid, fxciData.integrand[theta]);
      for (size_t i = 0; i < nrs; ++i) {
	const double &rs = rsGrid[i];
	if (rs <= rsMaxi) { fxcIntegrand[theta][i] = itp.eval(rs); }
      }
    }
  }
}

void ThermoPropBase::setAlpha(const VSInput &in) {
  // Initialize
  const size_t nrs = rsGrid.size();
  alpha.resize(nrs);
  // Set free parameter if passed in input
  const auto &fxciData = in.getFreeEnergyIntegrand();
  if (!fxciData.grid.empty()) {
    const double rsMaxi = fxciData.grid.back();
    for (size_t i = 0; i < nrs; ++i) {
      const double &rs = rsGrid[i];
      if (rs <= rsMaxi) { alpha[i] = fxciData.alpha[i]; }
    }
  }
}

void ThermoPropBase::setFxcIdxTargetStatePoint(const VSInput &in) {
  auto isTarget = [&](const double &rs) {
    return numUtil::equalTol(rs, in.getCoupling());
  };
  const auto it = std::find_if(rsGrid.begin(), rsGrid.end(), isTarget);
  if (it == rsGrid.end()) {
    MPIUtil::throwError("Failed to find the target state point in the free energy grid");
  }
  fxcIdxTargetStatePoint = std::distance(rsGrid.begin(), it);
}

void ThermoPropBase::setFxcIdxUnsolvedStatePoint() {
  const auto &fxciBegin = fxcIntegrand[Idx::THETA].begin();
  const auto &fxciEnd = fxcIntegrand[Idx::THETA].end();
  const auto &it = std::find(fxciBegin, fxciEnd, numUtil::Inf);
  fxcIdxUnsolvedStatePoint = std::distance(fxciBegin, it);
}


void ThermoPropBase::copyFreeEnergyIntegrand(const ThermoPropBase &other) {
  assert(other.rsGrid[1] - other.rsGrid[0] == rsGrid[1] - rsGrid[0]);
  const size_t nrs = rsGrid.size();
  const size_t nrsOther = other.rsGrid.size();
  for (const auto &theta : {Idx::THETA_DOWN, Idx::THETA, Idx::THETA_UP}) {
    const auto &fxciBegin = fxcIntegrand[theta].begin();
    const auto &fxciEnd = fxcIntegrand[theta].end();
    const auto &it = std::find(fxciBegin, fxciEnd, numUtil::Inf);
    size_t i = std::distance(fxciBegin, it);
    while (i < nrs && i < nrsOther) {
      fxcIntegrand[theta][i] = other.fxcIntegrand[theta][i];
      ++i;
    }
  }
  setFxcIdxUnsolvedStatePoint();
}

void ThermoPropBase::setAlpha(const double &alpha) {
  //structProp.setAlpha(alpha);
}

bool ThermoPropBase::isFreeEnergyIntegrandIncomplete() const {
  return fxcIdxUnsolvedStatePoint < fxcIdxTargetStatePoint - 1;
}

double ThermoPropBase::getFirstUnsolvedStatePoint() const {
  if (isFreeEnergyIntegrandIncomplete()) {
    return rsGrid[fxcIdxUnsolvedStatePoint + 1];
  } else {
    return numUtil::Inf;
  }
}

void ThermoPropBase::compute() {
  // structProp.compute();
  // const std::vector<double> fxciTmp = structProp.getFreeEnergyIntegrand();
  // const double alphaTmp = structProp.getAlpha();
  const std::vector<double> fxciTmp(3 * NPOINTS);
  const double alphaTmp = -1;
  const size_t &idx = fxcIdxTargetStatePoint;
  fxcIntegrand[THETA_DOWN][idx - 1] = fxciTmp[SIdx::RS_DOWN_THETA_DOWN];
  fxcIntegrand[THETA_DOWN][idx] = fxciTmp[SIdx::RS_THETA_DOWN];
  fxcIntegrand[THETA_DOWN][idx + 1] = fxciTmp[SIdx::RS_UP_THETA_DOWN];
  fxcIntegrand[THETA][idx - 1] = fxciTmp[SIdx::RS_DOWN_THETA];
  fxcIntegrand[THETA][idx] = fxciTmp[SIdx::RS_THETA];
  fxcIntegrand[THETA][idx + 1] = fxciTmp[SIdx::RS_UP_THETA];
  fxcIntegrand[THETA_UP][idx - 1] = fxciTmp[SIdx::RS_DOWN_THETA_UP];
  fxcIntegrand[THETA_UP][idx] = fxciTmp[SIdx::RS_THETA_UP];
  fxcIntegrand[THETA_UP][idx + 1] = fxciTmp[SIdx::RS_UP_THETA_UP];
  alpha[idx - 1] = alphaTmp;
  alpha[idx] = alphaTmp;
  alpha[idx + 1] = alphaTmp;
}

std::vector<double> ThermoPropBase::getSsf() {
  // if (!structProp.isComputed()) { structProp.compute(); }
  // return structProp.getCsr(getStructPropIdx()).getSsf();
  return std::vector<double>();
}

std::vector<double> ThermoPropBase::getSlfc() {
  // if (!structProp.isComputed()) { structProp.compute(); }
  // return structProp.getCsr(getStructPropIdx()).getSlfc();
  return std::vector<double>();
}

std::vector<double> ThermoPropBase::getFreeEnergyData() const {
  // const std::vector<double> rsVec = structProp.getCouplingParameters();
  // const std::vector<double> thetaVec = structProp.getDegeneracyParameters();
  // // Free energy
  // const double fxc = computeFreeEnergy(SIdx::RS_THETA, true);
  // // Free energy derivatives with respect to the coupling parameter
  // double fxcr;
  // double fxcrr;
  // {
  //   const double rs = rsVec[SIdx::RS_THETA];
  //   const double drs = rsVec[SIdx::RS_UP_THETA] - rsVec[SIdx::RS_THETA];
  //   const double f0 = computeFreeEnergy(SIdx::RS_UP_THETA, false);
  //   const double f1 = computeFreeEnergy(SIdx::RS_THETA, false);
  //   const double f2 = computeFreeEnergy(SIdx::RS_DOWN_THETA, false);
  //   fxcr = (f0 - f2) / (2.0 * drs * rs) - 2.0 * fxc;
  //   fxcrr = (f0 - 2.0 * f1 + f2) / (drs * drs) - 2.0 * fxc - 4.0 * fxcr;
  // }
  // // Free energy derivatives with respect to the degeneracy parameter
  // double fxct;
  // double fxctt;
  // {
  //   const double theta = thetaVec[SIdx::RS_THETA];
  //   const double theta2 = theta * theta;
  //   const double dt = thetaVec[SIdx::RS_THETA_UP] - thetaVec[SIdx::RS_THETA];
  //   const double f0 = computeFreeEnergy(SIdx::RS_THETA_UP, true);
  //   const double f1 = computeFreeEnergy(SIdx::RS_THETA_DOWN, true);
  //   fxct = theta * (f0 - f1) / (2.0 * dt);
  //   fxctt = theta2 * (f0 - 2.0 * fxc + f1) / (dt * dt);
  // }
  // // Free energy mixed derivatives
  // double fxcrt;
  // {
  //   const double t_rs = thetaVec[SIdx::RS_THETA] / rsVec[SIdx::RS_THETA];
  //   const double drs = rsVec[SIdx::RS_UP_THETA] - rsVec[SIdx::RS_THETA];
  //   const double dt = thetaVec[SIdx::RS_THETA_UP] - thetaVec[SIdx::RS_THETA];
  //   const double f0 = computeFreeEnergy(SIdx::RS_UP_THETA_UP, false);
  //   const double f1 = computeFreeEnergy(SIdx::RS_UP_THETA_DOWN, false);
  //   const double f2 = computeFreeEnergy(SIdx::RS_DOWN_THETA_UP, false);
  //   const double f3 = computeFreeEnergy(SIdx::RS_DOWN_THETA_DOWN, false);
  //   fxcrt = t_rs * (f0 - f1 - f2 + f3) / (4.0 * drs * dt) - 2.0 * fxct;
  // }
  // return std::vector<double>({fxc, fxcr, fxcrr, fxct, fxctt, fxcrt});
  return std::vector<double>();
}


std::vector<double> ThermoPropBase::getInternalEnergyData() const {
  // // Internal energy
  // const std::vector<double> uVec = structProp.getInternalEnergy();
  // const double u = uVec[SIdx::RS_THETA];
  // // Internal energy derivative with respect to the coupling parameter
  // double ur;
  // {
  //   const std::vector<double> rs = structProp.getCouplingParameters();
  //   const double drs = rs[SIdx::RS_UP_THETA] - rs[SIdx::RS_THETA];
  //   const std::vector<double> rsu = structProp.getFreeEnergyIntegrand();
  //   const double &u0 = rsu[SIdx::RS_UP_THETA];
  //   const double &u1 = rsu[SIdx::RS_DOWN_THETA];
  //   ur = (u0 - u1) / (2.0 * drs) - u;
  // }
  // // Internal energy derivative with respect to the degeneracy parameter
  // double ut;
  // {
  //   const std::vector<double> theta = structProp.getDegeneracyParameters();
  //   const double dt = theta[SIdx::RS_THETA_UP] - theta[SIdx::RS_THETA];
  //   const double u0 = uVec[SIdx::RS_THETA_UP];
  //   const double u1 = uVec[SIdx::RS_THETA_DOWN];
  //   ut = theta[SIdx::RS_THETA] * (u0 - u1) / (2.0 * dt);
  // }
  // return std::vector<double>({u, ur, ut});
  return std::vector<double>();
}


double ThermoPropBase::computeFreeEnergy(const ThermoPropBase::SIdx iStruct, const bool normalize) const {
  Idx iThermo;
  switch (iStruct) {
  case SIdx::RS_DOWN_THETA_DOWN:
  case SIdx::RS_THETA_DOWN:
  case SIdx::RS_UP_THETA_DOWN: iThermo = THETA_DOWN; break;
  case SIdx::RS_DOWN_THETA:
  case SIdx::RS_THETA:
  case SIdx::RS_UP_THETA: iThermo = THETA; break;
  case SIdx::RS_DOWN_THETA_UP:
  case SIdx::RS_THETA_UP:
  case SIdx::RS_UP_THETA_UP: iThermo = THETA_UP; break;
  default:
    assert(false);
    iThermo = THETA;
    break;
  }
  // const std::vector<double> &rs = structProp.getCouplingParameters();
  const std::vector<double> rs;
  return thermoUtil::computeFreeEnergy(rsGrid, fxcIntegrand[iThermo], rs[iStruct], normalize);
}

ThermoPropBase::SIdx ThermoPropBase::getStructPropIdx() {
  if (isZeroCoupling && isZeroDegeneracy) { return SIdx::RS_DOWN_THETA_DOWN; }
  if (!isZeroCoupling && isZeroDegeneracy) { return SIdx::RS_THETA_DOWN; }
  if (isZeroCoupling && !isZeroDegeneracy) { return SIdx::RS_DOWN_THETA; }
  return SIdx::RS_THETA;
}
