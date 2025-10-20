#include "vsbase.hpp"
#include "format.hpp"
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
  assert(thermoProp);
  return thermoProp->getFreeEnergyIntegrand();
}

const vector<double> &VSBase::getFreeEnergyGrid() const {
  assert(thermoProp);
  return thermoProp->getFreeEnergyGrid();
}

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
  assert(thermoProp);
  alpha = alphaTmp;
  thermoProp->setAlpha(alpha);
  const double alphaTheoretical = computeAlpha();
  return alpha - alphaTheoretical;
}

// -----------------------------------------------------------------
// ThermoPropBase class
// -----------------------------------------------------------------

ThermoPropBase::ThermoPropBase(const std::shared_ptr<const VSInput> &inPtr_)
    : inPtr(inPtr_) {
  // Check if we are solving for particular state points
  isZeroCoupling = (inRpa().getCoupling() == 0.0);
  isZeroDegeneracy = (inRpa().getDegeneracy() == 0.0);
  // Set variables related to the free energy calculation
  setRsGrid();
  setFxcIntegrand();
  setFxcIdxTargetStatePoint();
  setFxcIdxUnsolvedStatePoint();
}

void ThermoPropBase::setRsGrid() {
  const double &rs = inRpa().getCoupling();
  const double &drs = in().getCouplingResolution();
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

void ThermoPropBase::setFxcIntegrand() {
  const size_t nrs = rsGrid.size();
  // Initialize the free energy integrand
  fxcIntegrand.resize(NPOINTS);
  for (auto &f : fxcIntegrand) {
    f.resize(nrs);
    vecUtil::fill(f, numUtil::Inf);
  }
  // Fill the free energy integrand and the free parameter if passed in input
  const auto &fxciData = in().getFreeEnergyIntegrand();
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

void ThermoPropBase::setFxcIdxTargetStatePoint() {
  auto isTarget = [&](const double &rs) {
    return numUtil::equalTol(rs, inRpa().getCoupling());
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
  alpha = structProp->getAlpha();
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
}

const vector<double> &ThermoPropBase::getSsf() {
  assert(structProp);
  if (!structProp->isComputed()) { structProp->compute(); }
  std::cerr << "ThermoPropBase::getSsf() FIX INDEX ACCESSING " << std::endl;
  return structProp->getSsf(); // Fix correct index accessing
  // return structProp->getCsr(getStructPropIdx()).getSsf();
}

const Vector2D &ThermoPropBase::getLfc() {
  assert(structProp);
  if (!structProp->isComputed()) { structProp->compute(); }
  return structProp->getLfc(); // Fix correct index accessing
  // return structProp->getCsr(getStructPropIdx()).getLfc();
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

StructPropBase::StructPropBase(
    const std::shared_ptr<const IterationInput> &inPtr_)
    : inPtr(inPtr_),
      computed(false) {}

int StructPropBase::compute() {
  int status = csr->compute();
  println(formatUtil::format("Alpha = {:.5e}, Residual error "
                             "(structural properties) = {:.5e}",
                             csr->getAlpha(),
                             csr->getError()));
  computed = true;
  return status;
}

double StructPropBase::getAlpha() const { return csr->getAlpha(); }

void StructPropBase::setAlpha(const double &alpha) { csr->setAlpha(alpha); }

const vector<double> StructPropBase::getCouplingParameters() const {
  return csr->getAllCouplingParameters();
}

const vector<double> StructPropBase::getDegeneracyParameters() const {
  return csr->getAllDegeneracyParameters();
}

const vector<double> StructPropBase::getInternalEnergy() const {
  return csr->getAllInternalEnergies();
}

const vector<double> StructPropBase::getFreeEnergyIntegrand() const {
  return csr->getAllFreeEnergyIntegrands();
}

const vector<double> &StructPropBase::getSsf() const { return csr->getSsf(); }

const Vector2D &StructPropBase::getLfc() const { return csr->getLfc(); }

// -----------------------------------------------------------------
// CSRNew class
// -----------------------------------------------------------------

// Get the free parameter
double CSRNew::getAlpha() const {
  if (isManager) { return workers[StructIdx::RS_THETA]->getAlpha(); }
  return alpha;
}

void CSRNew::setAlpha(const double &alpha_) {
  if (isManager) {
    for (auto &worker : workers) {
      worker->setAlpha(alpha_);
    }
  } else {
    alpha = alpha_;
  }
}

std::vector<double> CSRNew::getAllCouplingParameters() const {
  vector<double> all;
  if (isManager) {
    for (auto &worker : workers) {
      all.push_back(worker->inRpa().getCoupling());
    }
  }
  return all;
}

std::vector<double> CSRNew::getAllDegeneracyParameters() const {
  vector<double> all;
  if (isManager) {
    for (auto &worker : workers) {
      all.push_back(worker->inRpa().getDegeneracy());
    }
  }
  return all;
}

std::vector<double> CSRNew::getAllInternalEnergies() const {
  vector<double> all;
  if (isManager) {
    for (auto &worker : workers) {
      const auto rs = worker->inRpa().getCoupling();
      const auto wvg = worker->getWvg();
      const auto ssf = worker->getSsf();
      const auto dim = worker->inRpa().getDimension();
      const double uint = thermoUtil::computeInternalEnergy(wvg, ssf, rs, dim);
      all.push_back(uint);
    }
  }
  return all;
}

std::vector<double> CSRNew::getAllFreeEnergyIntegrands() const {
  vector<double> all;
  if (isManager) {
    for (auto &worker : workers) {
      const auto wvg = worker->getWvg();
      const auto ssf = worker->getSsf();
      const auto dim = worker->inRpa().getDimension();
      const double uint = thermoUtil::computeInternalEnergy(wvg, ssf, 1.0, dim);
      all.push_back(uint);
    }
  }
  return all;
}

void CSRNew::computeLfcDerivative() {
  // Check that alpha has been set to a value that is not the default
  assert(alpha != DEFAULT_ALPHA);
  // Derivative contributions
  const double &rs = inRpa().getCoupling();
  const double &theta = inRpa().getDegeneracy();
  const double &dx = inRpa().getWaveVectorGridRes();
  const double &drs = inVS().getCouplingResolution();
  const double &dTheta = inVS().getDegeneracyResolution();
  const Vector2D &lfc = getLfc();
  const Vector2D &rsUp = *lfcRs.up;
  const Vector2D &rsDown = *lfcRs.down;
  const Vector2D &thetaUp = *lfcTheta.up;
  const Vector2D &thetaDown = *lfcTheta.down;
  const double a_drs = alpha * rs / (6.0 * drs);
  const double a_dx = alpha / (6.0 * dx);
  const double a_dt = alpha * theta / (3.0 * dTheta);
  const vector<double> &wvg = getWvg();
  const double nx = wvg.size();
  Vector2D &lfcd = lfcDerivative;
  assert(lfcd.size(0) == lfc.size(0) && lfcd.size(1) == lfc.size(1));
  for (size_t l = 0; l < lfc.size(1); ++l) {
    // Wave-vector derivative contribution
    lfcd(0, l) = a_dx * wvg[0] * getDerivative(lfc, l, 0, FORWARD);
    for (size_t i = 1; i < nx - 1; ++i) {
      lfcd(i, l) = a_dx * wvg[i] * getDerivative(lfc, l, i, CENTERED);
    }
    lfcd(nx - 1, l) =
        a_dx * wvg[nx - 1] * getDerivative(lfc, l, nx - 1, BACKWARD);
    // Coupling parameter contribution
    if (rs > 0.0) {
      for (size_t i = 0; i < nx; ++i) {
        lfcd(i, l) +=
            a_drs
            * getDerivative(lfc(i, l), rsUp(i, l), rsDown(i, l), lfcRs.type);
      }
    }
    // Degeneracy parameter contribution
    if (theta > 0.0) {
      for (size_t i = 0; i < nx; ++i) {
        lfcd(i, l) +=
            a_dt
            * getDerivative(
                lfc(i, l), thetaUp(i, l), thetaDown(i, l), lfcTheta.type);
      }
    }
  }
}

double CSRNew::getDerivative(const Vector2D &f,
                             const int &l,
                             const size_t &idx,
                             const Derivative &type) const {
  switch (type) {
  case BACKWARD:
    assert(idx >= 2);
    return CSRNew::getDerivative(f(idx, l), f(idx - 1, l), f(idx - 2, l), type);
    break;
  case CENTERED:
    assert(idx >= 1 && idx < f.size() - 1);
    return CSRNew::getDerivative(f(idx, l), f(idx + 1, l), f(idx - 1, l), type);
    break;
  case FORWARD:
    assert(idx < f.size() - 2);
    return CSRNew::getDerivative(f(idx, l), f(idx + 1, l), f(idx + 2, l), type);
    break;
  default:
    assert(false);
    return -1;
    break;
  }
}

double CSRNew::getDerivative(const double &f0,
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

void CSRNew::setupDerivativeData() {
  // if (!isManager) { return; }
  // for (size_t i = 0; i < workers.size(); ++i) {
  //   auto worker = workers[i];
  //   switch (i) {
  //   case RS_DOWN_THETA_DOWN:
  //     worker->setDrsData(*workers[RS_THETA_DOWN],
  //                        *workers[RS_UP_THETA_DOWN],
  //                        Derivative::FORWARD);
  //     break;
  //   case RS_THETA_DOWN:
  //     worker->setDrsData(*workers[RS_UP_THETA_DOWN],
  //                        *workers[RS_DOWN_THETA_DOWN],
  //                        Derivative::CENTERED);
  //     break;
  //   case RS_UP_THETA_DOWN:
  //     worker->setDrsData(*workers[RS_THETA_DOWN],
  //                        *workers[RS_DOWN_THETA_DOWN],
  //                        Derivative::BACKWARD);
  //     break;
  //   case RS_DOWN_THETA:
  //     worker->setDrsData(
  //         *workers[RS_THETA], *workers[RS_UP_THETA], Derivative::FORWARD);
  //     break;
  //   case RS_UP_THETA:
  //     worker->setDrsData(
  //         *workers[RS_THETA], *workers[RS_DOWN_THETA], Derivative::BACKWARD);
  //     break;
  //   case RS_DOWN_THETA_UP:
  //     worker->setDrsData(
  //         *workers[RS_THETA_UP], *workers[RS_UP_THETA_UP], Derivative::FORWARD);
  //     break;
  //   case RS_THETA_UP:
  //     worker->setDrsData(*workers[RS_UP_THETA_UP],
  //                        *workers[RS_DOWN_THETA_UP],
  //                        Derivative::CENTERED);
  //     break;
  //   case RS_UP_THETA_UP:
  //     worker->setDrsData(*workers[RS_THETA_UP],
  //                        *workers[RS_DOWN_THETA_UP],
  //                        Derivative::BACKWARD);
  //     break;
  //   }
  // }
  // for (size_t i = 0; i < workers.size(); ++i) {
  //   auto worker = workers[i];
  //   switch (i) {
  //   case RS_DOWN_THETA_DOWN:
  //     worker->setDThetaData(*workers[RS_DOWN_THETA],
  //                           *workers[RS_DOWN_THETA_UP],
  //                           Derivative::FORWARD);
  //     break;
  //   case RS_THETA_DOWN:
  //     worker->setDThetaData(
  //         *workers[RS_THETA], *workers[RS_THETA_UP], Derivative::FORWARD);
  //     break;
  //   case RS_UP_THETA_DOWN:
  //     worker->setDThetaData(
  //         *workers[RS_UP_THETA], *workers[RS_UP_THETA_UP], Derivative::FORWARD);
  //     break;
  //   case RS_DOWN_THETA:
  //     worker->setDThetaData(*workers[RS_DOWN_THETA_UP],
  //                           *workers[RS_DOWN_THETA_DOWN],
  //                           Derivative::CENTERED);
  //     break;
  //   case RS_UP_THETA:
  //     worker->setDThetaData(*workers[RS_UP_THETA_UP],
  //                           *workers[RS_UP_THETA_DOWN],
  //                           Derivative::CENTERED);
  //     break;
  //   case RS_DOWN_THETA_UP:
  //     worker->setDThetaData(*workers[RS_DOWN_THETA],
  //                           *workers[RS_DOWN_THETA_DOWN],
  //                           Derivative::BACKWARD);
  //     break;
  //   case RS_THETA_UP:
  //     worker->setDThetaData(
  //         *workers[RS_THETA], *workers[RS_THETA_DOWN], Derivative::BACKWARD);
  //     break;
  //   case RS_UP_THETA_UP:
  //     worker->setDThetaData(*workers[RS_UP_THETA],
  //                           *workers[RS_UP_THETA_DOWN],
  //                           Derivative::BACKWARD);
  //     break;
  //   }
  // }
  for (size_t i = 0; i < workers.size(); ++i) {
    auto &worker = *workers[i];
    switch (i) {
    case RS_DOWN_THETA_DOWN:
    case RS_DOWN_THETA:
    case RS_DOWN_THETA_UP:
      worker.setDrsData(*workers[i + 1], *workers[i + 2], Derivative::FORWARD);
      break;
    case RS_THETA_DOWN:
    case RS_THETA:
    case RS_THETA_UP:
      worker.setDrsData(*workers[i + 1], *workers[i - 1], Derivative::CENTERED);
      break;
    case RS_UP_THETA_DOWN:
    case RS_UP_THETA:
    case RS_UP_THETA_UP:
      worker.setDrsData(*workers[i - 1], *workers[i - 2], Derivative::BACKWARD);
      break;
    }
  }
  for (size_t i = 0; i < workers.size(); ++i) {
    auto &worker = *workers[i];
    switch (i) {
    case StructIdx::RS_DOWN_THETA_DOWN:
    case StructIdx::RS_THETA_DOWN:
    case StructIdx::RS_UP_THETA_DOWN:
      worker.setDThetaData(
          *workers[i + NRS], *workers[i + 2 * NRS], Derivative::FORWARD);
      break;
    case StructIdx::RS_DOWN_THETA:
    case StructIdx::RS_THETA:
    case StructIdx::RS_UP_THETA:
      worker.setDThetaData(
          *workers[i + NRS], *workers[i - NRS], Derivative::CENTERED);
      break;
    case StructIdx::RS_DOWN_THETA_UP:
    case StructIdx::RS_THETA_UP:
    case StructIdx::RS_UP_THETA_UP:
      worker.setDThetaData(
          *workers[i - NRS], *workers[i - 2 * NRS], Derivative::BACKWARD);
      break;
    }
  }
}

void CSRNew::setDrsData(CSRNew &up, CSRNew &down, const Derivative &dType) {
  lfcRs = DerivativeData{dType, &up.getLfc(), &down.getLfc()};
}

void CSRNew::setDThetaData(CSRNew &up, CSRNew &down, const Derivative &dType) {
  lfcTheta = DerivativeData{dType, &up.getLfc(), &down.getLfc()};
}

void CSRNew::init() {
  if (isManager) {
    for (auto &worker : workers) {
      worker->init();
    }
  } else {
    if (isInitialized) { return; }
    initStls();
    isInitialized = true;
  }
}

void CSRNew::initialGuess() {
  if (isManager) {
    for (auto &worker : workers) {
      worker->initialGuess();
    }
  } else {
    initialGuessStls();
  }
}

void CSRNew::computeSsf() {
  if (isManager) {
    for (auto &worker : workers) {
      worker->computeSsf();
    }
  } else {
    computeSsfStls();
  }
}

void CSRNew::computeLfc() {
  computeLfcPart1();
  computeLfcPart2();
  computeLfcPart3();
}

void CSRNew::computeLfcPart1() {
  if (isManager) {
    for (auto &worker : workers) {
      worker->computeLfcPart1();
    }
  } else {
    computeLfcStls();
    if (lfcDerivative.empty()) {
      lfcDerivative.resize(getLfc().size(0), getLfc().size(1));
    }
  }
}

void CSRNew::computeLfcPart2() {
  if (isManager) {
    for (auto &worker : workers) {
      worker->computeLfcPart2();
    }
  } else {
    computeLfcDerivative();
  }
}

void CSRNew::computeLfcPart3() {
  if (isManager) {
    for (auto &worker : workers) {
      worker->computeLfcPart3();
    }
  } else {
    getLfc().diff(lfcDerivative);
  }
}

void CSRNew::updateSolution() {
  if (isManager) {
    for (auto &worker : workers) {
      worker->updateSolution();
    }
  } else {
    updateSolutionStls();
  }
}

double CSRNew::computeError() const {
  if (isManager) { return workers[StructIdx::RS_THETA]->computeError(); }
  return computeErrorStls();
}