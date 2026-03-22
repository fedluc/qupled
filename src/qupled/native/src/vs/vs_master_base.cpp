#include "vs/vs_master_base.hpp"
#include "num_util.hpp"
#include "thermo_util.hpp"
#include <cassert>

using namespace std;
using namespace GridPoints;

// -----------------------------------------------------------------
// VSMasterBase
// -----------------------------------------------------------------

VSMasterBase::VSMasterBase(double drs_,
                           double dTheta_,
                           double dx_,
                           dimensionsUtil::Dimension dim_)
    : alpha(numUtil::Inf),
      lastError(numUtil::Inf),
      initDone(false),
      drs(drs_),
      dTheta(dTheta_),
      dx(dx_),
      dim(dim_) {}

void VSMasterBase::setupDerivativeData() {
  // rs derivative data (varies along rs axis, theta held fixed)
  for (const auto tOff : {GridPoint::Theta::DOWN,
                          GridPoint::Theta::CENTER,
                          GridPoint::Theta::UP}) {
    const size_t iDown = GridPoint{GridPoint::Rs::DOWN, tOff}.toIndex();
    const size_t iCenter = GridPoint{GridPoint::Rs::CENTER, tOff}.toIndex();
    const size_t iUp = GridPoint{GridPoint::Rs::UP, tOff}.toIndex();
    rsDerivData[iDown] = {DerivativeData::Type::FORWARD, iCenter, iUp};
    rsDerivData[iCenter] = {DerivativeData::Type::CENTERED, iUp, iDown};
    rsDerivData[iUp] = {DerivativeData::Type::BACKWARD, iCenter, iDown};
  }
  // theta derivative data (varies along theta axis, rs held fixed)
  for (const auto rOff :
       {GridPoint::Rs::DOWN, GridPoint::Rs::CENTER, GridPoint::Rs::UP}) {
    const size_t iDown = GridPoint{rOff, GridPoint::Theta::DOWN}.toIndex();
    const size_t iCenter = GridPoint{rOff, GridPoint::Theta::CENTER}.toIndex();
    const size_t iUp = GridPoint{rOff, GridPoint::Theta::UP}.toIndex();
    thetaDerivData[iDown] = {DerivativeData::Type::FORWARD, iCenter, iUp};
    thetaDerivData[iCenter] = {DerivativeData::Type::CENTERED, iUp, iDown};
    thetaDerivData[iUp] = {DerivativeData::Type::BACKWARD, iCenter, iDown};
  }
}

void VSMasterBase::computeSynchronizedLfc() {
  // Step 1: base LFC for every worker (must all finish before step 2)
  for (auto &w : workers) {
    w->computeBaseLfc();
  }
  // Resize derivative arrays on first call
  for (size_t i = 0; i < N; ++i) {
    const Vector2D &lfc = workers[i]->getLfc();
    if (lfcDerivatives[i].empty()) {
      lfcDerivatives[i].resize(lfc.size(0), lfc.size(1));
    }
  }
  // Step 2: compute the alpha-correction derivative for every worker
  computeLfcDerivatives();
  // Step 3: apply correction:  lfc -= lfcDerivative
  for (size_t i = 0; i < N; ++i) {
    workers[i]->applyLfcDiff(lfcDerivatives[i]);
  }
}

void VSMasterBase::computeLfcDerivatives() {
  assert(alpha != numUtil::Inf);

  for (size_t i = 0; i < N; ++i) {
    const double rs = rsValues[i];
    const double theta = thetaValues[i];
    const Vector2D &lfc = workers[i]->getLfc();
    const Vector2D &rsUpLfc = workers[rsDerivData[i].upIdx]->getLfc();
    const Vector2D &rsDownLfc = workers[rsDerivData[i].downIdx]->getLfc();
    const Vector2D &tUpLfc = workers[thetaDerivData[i].upIdx]->getLfc();
    const Vector2D &tDownLfc = workers[thetaDerivData[i].downIdx]->getLfc();
    const double a_dx = alpha / (6.0 * dx);
    const double a_drs = (rs > 0.0) ? alpha * rs / (6.0 * drs) : 0.0;
    const double a_dt = (theta > 0.0) ? alpha * theta / (3.0 * dTheta) : 0.0;
    const vector<double> &wvg = workers[i]->getWvg();
    const size_t nx = wvg.size();
    Vector2D &lfcd = lfcDerivatives[i];

    for (size_t l = 0; l < lfc.size(1); ++l) {
      // Wave-vector derivative
      lfcd(0, l) =
          a_dx * wvg[0] * derivative(lfc, l, 0, DerivativeData::Type::FORWARD);
      for (size_t k = 1; k < nx - 1; ++k) {
        lfcd(k, l) = a_dx * wvg[k]
                     * derivative(lfc, l, k, DerivativeData::Type::CENTERED);
      }
      lfcd(nx - 1, l) =
          a_dx * wvg[nx - 1]
          * derivative(lfc, l, nx - 1, DerivativeData::Type::BACKWARD);
      // Coupling parameter contribution
      for (size_t k = 0; k < nx; ++k) {
        lfcd(k, l) +=
            a_drs
            * derivative(
                lfc(k, l), rsUpLfc(k, l), rsDownLfc(k, l), rsDerivData[i].type);
      }
      // Degeneracy parameter contribution
      for (size_t k = 0; k < nx; ++k) {
        lfcd(k, l) += a_dt
                      * derivative(lfc(k, l),
                                   tUpLfc(k, l),
                                   tDownLfc(k, l),
                                   thetaDerivData[i].type);
      }
    }
  }
}

double VSMasterBase::derivative(const Vector2D &f,
                                int l,
                                size_t i,
                                DerivativeData::Type t) const {
  switch (t) {
  case DerivativeData::Type::BACKWARD:
    assert(i >= 2);
    return derivative(f(i, l), f(i - 1, l), f(i - 2, l), t);
  case DerivativeData::Type::CENTERED:
    assert(i >= 1 && i < f.size() - 1);
    return derivative(f(i, l), f(i + 1, l), f(i - 1, l), t);
  case DerivativeData::Type::FORWARD:
    assert(i < f.size() - 2);
    return derivative(f(i, l), f(i + 1, l), f(i + 2, l), t);
  default: assert(false); return -1.0;
  }
}

double VSMasterBase::derivative(double f0,
                                double f1,
                                double f2,
                                DerivativeData::Type t) const {
  switch (t) {
  case DerivativeData::Type::BACKWARD: return 3.0 * f0 - 4.0 * f1 + f2;
  case DerivativeData::Type::CENTERED: return f1 - f2;
  case DerivativeData::Type::FORWARD: return -(3.0 * f0 - 4.0 * f1 + f2);
  default: assert(false); return -1.0;
  }
}

// Iteration helpers

void VSMasterBase::masterInit() {
  if (initDone) return;
  for (auto &w : workers)
    w->init();
  initDone = true;
}

void VSMasterBase::masterComputeLfc() { computeSynchronizedLfc(); }

void VSMasterBase::masterComputeSsf() {
  for (auto &w : workers)
    w->computeSsf();
}

double VSMasterBase::masterComputeError() const {
  lastError = workers[CENTER.toIndex()]->computeError();
  return lastError;
}

void VSMasterBase::masterUpdateSolution() {
  for (auto &w : workers)
    w->updateSolution();
}

void VSMasterBase::masterInitialGuess() {
  for (auto &w : workers)
    w->initialGuess();
}

// Getters

const std::vector<double> &VSMasterBase::getSsf(GridPoint p) const {
  return workers[p.toIndex()]->getSsf();
}

const Vector2D &VSMasterBase::getLfc(GridPoint p) const {
  return workers[p.toIndex()]->getLfc();
}

const std::vector<double> &VSMasterBase::getWvg(GridPoint p) const {
  return workers[p.toIndex()]->getWvg();
}

double VSMasterBase::getCoupling(GridPoint p) const {
  return rsValues[p.toIndex()];
}

double VSMasterBase::getDegeneracy(GridPoint p) const {
  return thetaValues[p.toIndex()];
}

double VSMasterBase::getUInt(GridPoint p) const {
  const size_t i = p.toIndex();
  return thermoUtil::computeInternalEnergy(
      workers[i]->getWvg(), workers[i]->getSsf(), rsValues[i], dim);
}

double VSMasterBase::getFxcIntegrandValue(GridPoint p) const {
  const size_t i = p.toIndex();
  return thermoUtil::computeInternalEnergy(
      workers[i]->getWvg(), workers[i]->getSsf(), 1.0, dim);
}

const VSWorkerBase &VSMasterBase::getWorkerAt(GridPoint p) const {
  return *workers[p.toIndex()];
}
