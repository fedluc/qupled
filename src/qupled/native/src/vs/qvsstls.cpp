#include "vs/qvsstls.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;
using namespace MPIUtil;
using namespace GridPoints;
using ItgType = Integrator1D::Type;
using Itg2DParam = Integrator2D::Param;

// -----------------------------------------------------------------
// VSQstlsWorker
// -----------------------------------------------------------------

VSQstlsWorker::VSQstlsWorker(const std::shared_ptr<const QVSStlsInput> &in,
                             GridPoint p)
    : Qstls(in, true) {
  const string &theory = in->getTheory();
  switch (p.theta) {
  case GridPoint::Theta::DOWN:
    adrFixedDatabaseName = formatUtil::format("{}_THETA_DOWN", theory);
    break;
  case GridPoint::Theta::CENTER:
    adrFixedDatabaseName = formatUtil::format("{}_THETA", theory);
    break;
  case GridPoint::Theta::UP:
    adrFixedDatabaseName = formatUtil::format("{}_THETA_UP", theory);
    break;
  }
}

double VSQstlsWorker::computeQAdder(const std::shared_ptr<Integrator2D> &itg2D,
                                    const std::vector<double> &itgGrid) const {
  const auto ssfItp = make_shared<Interpolator1D>(wvg, ssf);
  const QAdder q = QAdder::quantum(inPtr->getDegeneracy(),
                                   mu,
                                   wvg.front(),
                                   wvg.back(),
                                   itgGrid,
                                   itg,
                                   itg2D,
                                   ssfItp);
  return q.get();
}

// -----------------------------------------------------------------
// VSQstlsManager
// -----------------------------------------------------------------

VSQstlsManager::VSQstlsManager(const std::shared_ptr<const QVSStlsInput> &in)
    : VSManager(in->getCouplingResolution(),
                in->getDegeneracyResolution(),
                in->getWaveVectorGridRes(),
                in->getDimension()),
      Qstls(in, false),
      itg2D(make_shared<Integrator2D>(ItgType::DEFAULT, in->getIntError())) {
  const double drs_ = in->getCouplingResolution();
  const double dTheta_ = in->getDegeneracyResolution();
  const double rs0 = std::max(in->getCoupling(), drs_);
  const double theta0 = std::max(in->getDegeneracy(), dTheta_);
  for (const auto tOff : {GridPoint::Theta::DOWN,
                          GridPoint::Theta::CENTER,
                          GridPoint::Theta::UP}) {
    const double thetaTmp = theta0 + static_cast<int>(tOff) * dTheta_;
    for (const auto rOff :
         {GridPoint::Rs::DOWN, GridPoint::Rs::CENTER, GridPoint::Rs::UP}) {
      const double rsTmp = rs0 + static_cast<int>(rOff) * drs_;
      auto inTmp = std::make_shared<QVSStlsInput>(*in);
      inTmp->setCoupling(rsTmp);
      inTmp->setDegeneracy(thetaTmp);
      const GridPoint gp{rOff, tOff};
      const size_t idx = gp.toIndex();
      rsValues[idx] = rsTmp;
      thetaValues[idx] = thetaTmp;
      workers[idx] = std::make_unique<VSQstlsWorker>(inTmp, gp);
    }
  }
  setupDerivativeData();
  // Setup integration grid
  const bool segregatedItg = in->getInt2DScheme() == "segregated";
  if (segregatedItg) { itgGrid = VSManager::getWvg(GridPoints::CENTER); }
}

double VSQstlsManager::computeQRaw(GridPoint p) const {
  const auto &w = dynamic_cast<const VSQstlsWorker &>(getWorkerAt(p));
  return w.computeQAdder(itg2D, itgGrid);
}

// -----------------------------------------------------------------
// QVSStls
// -----------------------------------------------------------------

QVSStls::QVSStls(const std::shared_ptr<const QVSStlsInput> &in)
    : VSBase(),
      inPtr(in),
      grid_(in) {
  if (in->getDegeneracy() == 0.0) {
    throwError(
        "Ground state calculations are not implemented for this scheme.");
  }
  setRsGrid();
  setFxcIntegrand();
}

const VSInput &QVSStls::in() const {
  return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
}

const Input &QVSStls::inScheme() const { return *inPtr; }

double QVSStls::computeQRaw(GridPoint p) const {
  return grid_.computeQRaw(p);
}
