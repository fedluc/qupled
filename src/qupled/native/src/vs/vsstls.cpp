#include "vs/vsstls.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "thermo_util.hpp"

using namespace std;
using namespace MPIUtil;
using namespace GridPoints;

// -----------------------------------------------------------------
// VSStlsMaster
// -----------------------------------------------------------------

VSStlsMaster::VSStlsMaster(const std::shared_ptr<const VSStlsInput> &in)
    : VSMasterBase(in->getCouplingResolution(),
                   in->getDegeneracyResolution(),
                   in->getWaveVectorGridRes(),
                   in->getDimension()),
      Stls(in, false) {
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
      auto inTmp = std::make_shared<VSStlsInput>(*in);
      inTmp->setCoupling(rsTmp);
      inTmp->setDegeneracy(thetaTmp);
      const GridPoint gp{rOff, tOff};
      const size_t idx = gp.toIndex();
      rsValues[idx] = rsTmp;
      thetaValues[idx] = thetaTmp;
      workers[idx] = std::make_unique<VSStlsWorker>(inTmp, false, gp);
    }
  }
  setupDerivativeData();
}

// -----------------------------------------------------------------
// VSStls
// -----------------------------------------------------------------

VSStls::VSStls(const std::shared_ptr<const VSStlsInput> &in)
    : VSBase(),
      inPtr(in),
      grid(in) {
  if (in->getDimension() == dimensionsUtil::Dimension::D2) {
    throwError("2D calculations are not implemented for this scheme.");
  }
  setRsGrid();
  setFxcIntegrand();
}

const VSInput &VSStls::in() const {
  return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
}

const Input &VSStls::inScheme() const { return *inPtr; }

int VSStls::runGrid() {
  grid.setAlpha(alpha);
  int status = grid.compute();
  println(formatUtil::format("Alpha = {:.5e}, Residual error "
                             "(structural properties) = {:.5e}",
                             grid.getAlpha(),
                             grid.VSMasterBase::getError()));
  return status;
}

double VSStls::getCoupling(GridPoint p) const { return grid.getCoupling(p); }

double VSStls::getDegeneracy(GridPoint p) const {
  return grid.getDegeneracy(p);
}

double VSStls::getFxcIntegrandValue(GridPoint p) const {
  return grid.getFxcIntegrandValue(p);
}

double VSStls::computeQRaw(GridPoint p) const {
  return QAdder::classical(
             grid.VSMasterBase::getWvg(p), grid.VSMasterBase::getSsf(p), inPtr)
      .get();
}

const std::vector<double> &VSStls::getSsf() const {
  return grid.getWorkerAt(GridPoints::CENTER).getSsf();
}

const Vector2D &VSStls::getLfc() const {
  return grid.getWorkerAt(GridPoints::CENTER).getLfc();
}

const vector<double> &VSStls::getWvg() const {
  return grid.getWorkerAt(GridPoints::CENTER).getWvg();
}

const Vector2D &VSStls::getIdr() const {
  return dynamic_cast<const Stls &>(grid.getWorkerAt(GridPoints::CENTER))
      .getIdr();
}

vector<double> VSStls::getSdr() const {
  return dynamic_cast<const Stls &>(grid.getWorkerAt(GridPoints::CENTER))
      .getSdr();
}

double VSStls::getUInt() const {
  return dynamic_cast<const Stls &>(grid.getWorkerAt(GridPoints::CENTER))
      .getUInt();
}

double VSStls::getError() const { return grid.VSMasterBase::getError(); }
