#include "vsbase.hpp"
#include "input.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"
#include <fmt/core.h>

using namespace std;

// -----------------------------------------------------------------
// VSBase class
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

void VSBase::init() {
  initScheme();
  initFreeEnergyIntegrand();
}

const vector<vector<double>> &VSBase::getFreeEnergyIntegrand() const {
  assert(thermoProp);
  return thermoProp->getFreeEnergyIntegrand();
}

const vector<double> &VSBase::getFreeEnergyGrid() const {
  assert(thermoProp);
  return thermoProp->getFreeEnergyGrid();
}

const vector<double> &VSBase::getAlpha() const {
  assert(thermoProp);
  return thermoProp->getAlpha();
}

void VSBase::doIterations() {
  auto func = [this](const double &alphaTmp) -> double {
    return alphaDifference(alphaTmp);
  };
  SecantSolver rsol(in.getErrMinAlpha(), in.getNIterAlpha());
  rsol.solve(func, in.getAlphaGuess());
  alpha = rsol.getSolution();
  println(fmt::format("Free parameter = {:.5f}", alpha));
  updateSolution();
}

double VSBase::alphaDifference(const double &alphaTmp) {
  assert(thermoProp);
  alpha = alphaTmp;
  thermoProp->setAlpha(alpha);
  const double alphaTheoretical = computeAlpha();
  return alpha - alphaTheoretical;
}

// -----------------------------------------------------------------
// ThermoPropBase class
// -----------------------------------------------------------------

ThermoPropBase::ThermoPropBase(const VSInput &inVS, const RpaInput &inRpa) {
  // Check if we are solving for particular state points
  isZeroCoupling = (inRpa.getCoupling() == 0.0);
  isZeroDegeneracy = (inRpa.getDegeneracy() == 0.0);
  // Set variables related to the free energy calculation
  setRsGrid(inVS, inRpa);
  setFxcIntegrand(inVS);
  setAlpha(inVS);
  setFxcIdxTargetStatePoint(inRpa);
  setFxcIdxUnsolvedStatePoint();
}

void ThermoPropBase::setRsGrid(const VSInput &inVS, const RpaInput &inRpa) {
  const double &rs = inRpa.getCoupling();
  const double &drs = inVS.getCouplingResolution();
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
      if (itp.isValid()) {
        for (size_t i = 0; i < nrs; ++i) {
          const double &rs = rsGrid[i];
          if (rs <= rsMaxi) { fxcIntegrand[theta][i] = itp.eval(rs); }
        }
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

void ThermoPropBase::setFxcIdxTargetStatePoint(const RpaInput &in) {
  auto isTarget = [&](const double &rs) {
    return numUtil::equalTol(rs, in.getCoupling());
  };
  const auto it = find_if(rsGrid.begin(), rsGrid.end(), isTarget);
  if (it == rsGrid.end()) {
    MPIUtil::throwError(
        "Failed to find the target state point in the free energy grid");
  }
  fxcIdxTargetStatePoint = distance(rsGrid.begin(), it);
}

void ThermoPropBase::setFxcIdxUnsolvedStatePoint() {
  const auto &fxciBegin = fxcIntegrand[Idx::THETA].begin();
  const auto &fxciEnd = fxcIntegrand[Idx::THETA].end();
  const auto &it = find(fxciBegin, fxciEnd, numUtil::Inf);
  fxcIdxUnsolvedStatePoint = distance(fxciBegin, it);
}

void ThermoPropBase::copyFreeEnergyIntegrand(const ThermoPropBase &other) {
  assert(other.rsGrid[1] - other.rsGrid[0] == rsGrid[1] - rsGrid[0]);
  const size_t nrs = rsGrid.size();
  const size_t nrsOther = other.rsGrid.size();
  for (const auto &theta : {Idx::THETA_DOWN, Idx::THETA, Idx::THETA_UP}) {
    const auto &fxciBegin = fxcIntegrand[theta].begin();
    const auto &fxciEnd = fxcIntegrand[theta].end();
    const auto &it = find(fxciBegin, fxciEnd, numUtil::Inf);
    size_t i = distance(fxciBegin, it);
    while (i < nrs && i < nrsOther) {
      fxcIntegrand[theta][i] = other.fxcIntegrand[theta][i];
      ++i;
    }
  }
  setFxcIdxUnsolvedStatePoint();
}

void ThermoPropBase::setAlpha(const double &alpha) {
  assert(structProp);
  structProp->setAlpha(alpha);
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
  assert(structProp);
  structProp->compute();
  const vector<double> fxciTmp = structProp->getFreeEnergyIntegrand();
  const double alphaTmp = structProp->getAlpha();
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

const vector<double> &ThermoPropBase::getSsf() {
  assert(structProp);
  if (!structProp->isComputed()) { structProp->compute(); }
  return structProp->getCsr(getStructPropIdx()).getSsf();
}

const vector<double> &ThermoPropBase::getSlfc() {
  assert(structProp);
  if (!structProp->isComputed()) { structProp->compute(); }
  return structProp->getCsr(getStructPropIdx()).getSlfc();
}

vector<double> ThermoPropBase::getFreeEnergyData() const {
  assert(structProp);
  const vector<double> rsVec = structProp->getCouplingParameters();
  const vector<double> thetaVec = structProp->getDegeneracyParameters();
  // Free energy
  const double fxc = computeFreeEnergy(SIdx::RS_THETA, true);
  // Free energy derivatives with respect to the coupling parameter
  double fxcr;
  double fxcrr;
  {
    const double rs = rsVec[SIdx::RS_THETA];
    const double drs = rsVec[SIdx::RS_UP_THETA] - rsVec[SIdx::RS_THETA];
    const double f0 = computeFreeEnergy(SIdx::RS_UP_THETA, false);
    const double f1 = computeFreeEnergy(SIdx::RS_THETA, false);
    const double f2 = computeFreeEnergy(SIdx::RS_DOWN_THETA, false);
    fxcr = (f0 - f2) / (2.0 * drs * rs) - 2.0 * fxc;
    fxcrr = (f0 - 2.0 * f1 + f2) / (drs * drs) - 2.0 * fxc - 4.0 * fxcr;
  }
  // Free energy derivatives with respect to the degeneracy parameter
  double fxct;
  double fxctt;
  {
    const double theta = thetaVec[SIdx::RS_THETA];
    const double theta2 = theta * theta;
    const double dt = thetaVec[SIdx::RS_THETA_UP] - thetaVec[SIdx::RS_THETA];
    const double f0 = computeFreeEnergy(SIdx::RS_THETA_UP, true);
    const double f1 = computeFreeEnergy(SIdx::RS_THETA_DOWN, true);
    fxct = theta * (f0 - f1) / (2.0 * dt);
    fxctt = theta2 * (f0 - 2.0 * fxc + f1) / (dt * dt);
  }
  // Free energy mixed derivatives
  double fxcrt;
  {
    const double t_rs = thetaVec[SIdx::RS_THETA] / rsVec[SIdx::RS_THETA];
    const double drs = rsVec[SIdx::RS_UP_THETA] - rsVec[SIdx::RS_THETA];
    const double dt = thetaVec[SIdx::RS_THETA_UP] - thetaVec[SIdx::RS_THETA];
    const double f0 = computeFreeEnergy(SIdx::RS_UP_THETA_UP, false);
    const double f1 = computeFreeEnergy(SIdx::RS_UP_THETA_DOWN, false);
    const double f2 = computeFreeEnergy(SIdx::RS_DOWN_THETA_UP, false);
    const double f3 = computeFreeEnergy(SIdx::RS_DOWN_THETA_DOWN, false);
    fxcrt = t_rs * (f0 - f1 - f2 + f3) / (4.0 * drs * dt) - 2.0 * fxct;
  }
  return vector<double>({fxc, fxcr, fxcrr, fxct, fxctt, fxcrt});
}

vector<double> ThermoPropBase::getInternalEnergyData() const {
  assert(structProp);
  // Internal energy
  const vector<double> uVec = structProp->getInternalEnergy();
  const double u = uVec[SIdx::RS_THETA];
  // Internal energy derivative with respect to the coupling parameter
  double ur;
  {
    const vector<double> rs = structProp->getCouplingParameters();
    const double drs = rs[SIdx::RS_UP_THETA] - rs[SIdx::RS_THETA];
    const vector<double> rsu = structProp->getFreeEnergyIntegrand();
    const double &u0 = rsu[SIdx::RS_UP_THETA];
    const double &u1 = rsu[SIdx::RS_DOWN_THETA];
    ur = (u0 - u1) / (2.0 * drs) - u;
  }
  // Internal energy derivative with respect to the degeneracy parameter
  double ut;
  {
    const vector<double> theta = structProp->getDegeneracyParameters();
    const double dt = theta[SIdx::RS_THETA_UP] - theta[SIdx::RS_THETA];
    const double u0 = uVec[SIdx::RS_THETA_UP];
    const double u1 = uVec[SIdx::RS_THETA_DOWN];
    ut = theta[SIdx::RS_THETA] * (u0 - u1) / (2.0 * dt);
  }
  return vector<double>({u, ur, ut});
}

double ThermoPropBase::computeFreeEnergy(const ThermoPropBase::SIdx iStruct,
                                         const bool normalize) const {
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
  assert(structProp);
  const vector<double> &rs = structProp->getCouplingParameters();
  return thermoUtil::computeFreeEnergy(
      rsGrid, fxcIntegrand[iThermo], rs[iStruct], normalize);
}

ThermoPropBase::SIdx ThermoPropBase::getStructPropIdx() {
  if (isZeroCoupling && isZeroDegeneracy) { return SIdx::RS_DOWN_THETA_DOWN; }
  if (!isZeroCoupling && isZeroDegeneracy) { return SIdx::RS_THETA_DOWN; }
  if (isZeroCoupling && !isZeroDegeneracy) { return SIdx::RS_DOWN_THETA; }
  return SIdx::RS_THETA;
}

// -----------------------------------------------------------------
// StructPropBase class
// -----------------------------------------------------------------

StructPropBase::StructPropBase()
    : csrIsInitialized(false),
      computed(false),
      outVector(NPOINTS) {}

void StructPropBase::setupCSRDependencies() {
  for (size_t i = 0; i < csr.size(); ++i) {
    switch (i) {
    case RS_DOWN_THETA_DOWN:
    case RS_DOWN_THETA:
    case RS_DOWN_THETA_UP:
      csr[i]->setDrsData(*csr[i + 1], *csr[i + 2], CSR::Derivative::FORWARD);
      break;
    case RS_THETA_DOWN:
    case RS_THETA:
    case RS_THETA_UP:
      csr[i]->setDrsData(*csr[i + 1], *csr[i - 1], CSR::Derivative::CENTERED);
      break;
    case RS_UP_THETA_DOWN:
    case RS_UP_THETA:
    case RS_UP_THETA_UP:
      csr[i]->setDrsData(*csr[i - 1], *csr[i - 2], CSR::Derivative::BACKWARD);
      break;
    }
  }
  for (size_t i = 0; i < csr.size(); ++i) {
    switch (i) {
    case RS_DOWN_THETA_DOWN:
    case RS_THETA_DOWN:
    case RS_UP_THETA_DOWN:
      csr[i]->setDThetaData(
          *csr[i + NRS], *csr[i + 2 * NRS], CSR::Derivative::FORWARD);
      break;
    case RS_DOWN_THETA:
    case RS_THETA:
    case RS_UP_THETA:
      csr[i]->setDThetaData(
          *csr[i + NRS], *csr[i - NRS], CSR::Derivative::CENTERED);
      break;
    case RS_DOWN_THETA_UP:
    case RS_THETA_UP:
    case RS_UP_THETA_UP:
      csr[i]->setDThetaData(
          *csr[i - NRS], *csr[i - 2 * NRS], CSR::Derivative::BACKWARD);
      break;
    }
  }
}

int StructPropBase::compute() {
  try {
    if (!csrIsInitialized) {
      assert(!csr.empty());
      for (auto &c : csr) {
        c->init();
      }
      csrIsInitialized = true;
    }
    doIterations();
    computed = true;
    return 0;
  } catch (const runtime_error &err) {
    cerr << err.what() << endl;
    return 1;
  }
}

void StructPropBase::setAlpha(const double &alpha) {
  for (auto &c : csr) {
    c->setAlpha(alpha);
  }
}

const CSR &StructPropBase::getCsr(const Idx &idx) const {
  assert(!csr.empty());
  return *csr[idx];
}

double StructPropBase::getAlpha() const {
  return getCsr(Idx::RS_THETA).getAlpha();
}

const vector<double> &
StructPropBase::getBase(function<double(const CSR &)> f) const {
  for (size_t i = 0; i < csr.size(); ++i) {
    outVector[i] = f(*csr[i]);
  }
  return outVector;
}

const vector<double> &StructPropBase::getCouplingParameters() const {
  return getBase([&](const CSR &c) { return c.getCoupling(); });
}

const vector<double> &StructPropBase::getDegeneracyParameters() const {
  return getBase([&](const CSR &c) { return c.getDegeneracy(); });
}

const vector<double> &StructPropBase::getInternalEnergy() const {
  return getBase([&](const CSR &c) { return c.getInternalEnergy(); });
}

const vector<double> &StructPropBase::getFreeEnergyIntegrand() const {
  return getBase([&](const CSR &c) { return c.getFreeEnergyIntegrand(); });
}

// -----------------------------------------------------------------
// CSR class
// -----------------------------------------------------------------

void CSR::setDrsData(CSR &csrRsUp, CSR &csrRsDown, const Derivative &dTypeRs) {
  lfcRs = DerivativeData{dTypeRs, csrRsUp.lfc, csrRsDown.lfc};
}

void CSR::setDThetaData(CSR &csrThetaUp,
                        CSR &csrThetaDown,
                        const Derivative &dTypeTheta) {
  lfcTheta = DerivativeData{dTypeTheta, csrThetaUp.lfc, csrThetaDown.lfc};
}

double CSR::getInternalEnergy() const {
  const double rs = inRpa.getCoupling();
  return thermoUtil::computeInternalEnergy(getWvg(), getSsf(), rs);
}

double CSR::getFreeEnergyIntegrand() const {
  return thermoUtil::computeInternalEnergy(getWvg(), getSsf(), 1.0);
}

Vector2D CSR::getDerivativeContribution() const {
  // Check that alpha has been set to a value that is not the default
  assert(alpha != DEFAULT_ALPHA);
  // Derivative contributions
  const double &rs = inRpa.getCoupling();
  const double &theta = inRpa.getDegeneracy();
  const double &dx = inRpa.getWaveVectorGridRes();
  const double &drs = inVS.getCouplingResolution();
  const double &dTheta = inVS.getDegeneracyResolution();
  const Vector2D &lfcData = *lfc;
  const Vector2D &rsUp = *lfcRs.up;
  const Vector2D &rsDown = *lfcRs.down;
  const Vector2D &thetaUp = *lfcTheta.up;
  const Vector2D &thetaDown = *lfcTheta.down;
  const double a_drs = alpha * rs / (6.0 * drs);
  const double a_dx = alpha / (6.0 * dx);
  const double a_dt = alpha * theta / (3.0 * dTheta);
  const vector<double> &wvg = getWvg();
  const double nx = wvg.size();
  Vector2D out(lfc->size(0), lfc->size(1));
  for (size_t l = 0; l < lfc->size(1); ++l) {
    // Wave-vector derivative contribution
    out(0, l) = a_dx * wvg[0] * getDerivative(lfc, l, 0, FORWARD);
    for (size_t i = 1; i < nx - 1; ++i) {
      out(i, l) = a_dx * wvg[i] * getDerivative(lfc, l, i, CENTERED);
    }
    out(nx - 1, l) =
        a_dx * wvg[nx - 1] * getDerivative(lfc, l, nx - 1, BACKWARD);
    // Coupling parameter contribution
    if (rs > 0.0) {
      for (size_t i = 0; i < nx; ++i) {
        out(i, l) += a_drs
                     * CSR::getDerivative(
                         lfcData(i, l), rsUp(i, l), rsDown(i, l), lfcRs.type);
      }
    }
    // Degeneracy parameter contribution
    if (theta > 0.0) {
      for (size_t i = 0; i < nx; ++i) {
        out(i, l) +=
            a_dt
            * CSR::getDerivative(
                lfcData(i, l), thetaUp(i, l), thetaDown(i, l), lfcTheta.type);
      }
    }
  }
  return out;
}

double CSR::getDerivative(const shared_ptr<Vector2D> &f,
                          const int &l,
                          const size_t &idx,
                          const Derivative &type) const {
  const Vector2D &fData = *f;
  switch (type) {
  case BACKWARD:
    assert(idx >= 2);
    return CSR::getDerivative(
        fData(idx, l), fData(idx - 1, l), fData(idx - 2, l), type);
    break;
  case CENTERED:
    assert(idx >= 1 && idx < fData.size() - 1);
    return CSR::getDerivative(
        fData(idx, l), fData(idx + 1, l), fData(idx - 1, l), type);
    break;
  case FORWARD:
    assert(idx < fData.size() - 2);
    return CSR::getDerivative(
        fData(idx, l), fData(idx + 1, l), fData(idx + 2, l), type);
    break;
  default:
    assert(false);
    return -1;
    break;
  }
}

double CSR::getDerivative(const double &f0,
                          const double &f1,
                          const double &f2,
                          const Derivative &type) const {
  switch (type) {
  case BACKWARD: return 3.0 * f0 - 4.0 * f1 + f2; break;
  case CENTERED: return f1 - f2; break;
  case FORWARD: return -getDerivative(f0, f1, f2, BACKWARD); break;
  default:
    assert(false);
    return -1;
    break;
  }
}
