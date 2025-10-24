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

void ThermoPropBase::setAlpha(const double &alpha) {
  assert(structProp);
  structProp->setAlpha(alpha);
}

void ThermoPropBase::compute() {
  assert(structProp);
  structProp->compute();
  println(formatUtil::format("Alpha = {:.5e}, Residual error "
                             "(structural properties) = {:.5e}",
                             structProp->getAlpha(),
                             structProp->getError()));
  const vector<double> fxciTmp = structProp->getAllFreeEnergyIntegrands();
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
  return structProp->getSsf(getStructPropIdx());
}

const Vector2D &ThermoPropBase::getLfc() {
  assert(structProp);
  return structProp->getLfc(getStructPropIdx());
}

vector<double> ThermoPropBase::getFreeEnergyData() const {
  assert(structProp);
  const vector<double> rsVec = structProp->getAllCouplingParameters();
  const vector<double> thetaVec = structProp->getAllDegeneracyParameters();
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
  const vector<double> uVec = structProp->getAllInternalEnergies();
  const double u = uVec[SIdx::RS_THETA];
  // Internal energy derivative with respect to the coupling parameter
  double ur;
  {
    const vector<double> rs = structProp->getAllCouplingParameters();
    const double drs = rs[SIdx::RS_UP_THETA] - rs[SIdx::RS_THETA];
    const vector<double> rsu = structProp->getAllFreeEnergyIntegrands();
    const double &u0 = rsu[SIdx::RS_UP_THETA];
    const double &u1 = rsu[SIdx::RS_DOWN_THETA];
    ur = (u0 - u1) / (2.0 * drs) - u;
  }
  // Internal energy derivative with respect to the degeneracy parameter
  double ut;
  {
    const vector<double> theta = structProp->getAllDegeneracyParameters();
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
  const vector<double> &rs = structProp->getAllCouplingParameters();
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
// CSRNew class
// -----------------------------------------------------------------

void CSRNew::forEachWorker(const function<void(CSRNew &)> &func) {
  if (isManager) {
    for (auto &worker : workers)
      worker->forEachWorker(func);
  } else {
    func(*this);
  }
}

vector<double>
CSRNew::collectFromWorkers(const function<double(const CSRNew &)> &func) const {
  std::vector<double> results;
  if (isManager) {
    for (auto &worker : workers) {
      auto sub = worker->collectFromWorkers(func);
      results.insert(results.end(), sub.begin(), sub.end());
    }
  } else {
    results.push_back(func(*this));
  }
  return results;
}

double CSRNew::getAlpha() const {
  if (isManager) { return workers[StructIdx::RS_THETA]->getAlpha(); }
  return alpha;
}

const std::vector<double> &CSRNew::getSsf(const size_t &idx) const {
  if (isManager) { return workers[idx]->getSsf(); }
  return getSsf();
}

const Vector2D &CSRNew::getLfc(const size_t &idx) const {
  if (isManager) { return workers[idx]->getLfc(); }
  return getLfc();
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
  auto func = [](const CSRNew &self) { return self.inRpa().getCoupling(); };
  return collectFromWorkers(func);
}

std::vector<double> CSRNew::getAllDegeneracyParameters() const {
  auto func = [](const CSRNew &self) { return self.inRpa().getDegeneracy(); };
  return collectFromWorkers(func);
}

std::vector<double> CSRNew::getAllInternalEnergies() const {
  auto func = [](const CSRNew &self) {
    const auto rs = self.inRpa().getCoupling();
    const auto wvg = self.getWvg();
    const auto ssf = self.getSsf();
    const auto dim = self.inRpa().getDimension();
    return thermoUtil::computeInternalEnergy(wvg, ssf, rs, dim);
  };
  return collectFromWorkers(func);
}

std::vector<double> CSRNew::getAllFreeEnergyIntegrands() const {
  auto func = [](const CSRNew &self) {
    const auto wvg = self.getWvg();
    const auto ssf = self.getSsf();
    const auto dim = self.inRpa().getDimension();
    return thermoUtil::computeInternalEnergy(wvg, ssf, 1.0, dim);
  };
  return collectFromWorkers(func);
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