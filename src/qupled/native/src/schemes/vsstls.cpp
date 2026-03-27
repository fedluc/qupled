#include "schemes/vsstls.hpp"
#include "util/format.hpp"
#include "schemes/input.hpp"
#include "util/mpi_util.hpp"
#include "thermo/thermo_util.hpp"

using namespace std;
using namespace MPIUtil;
using namespace GridPoints;

// -----------------------------------------------------------------
// VSStlsManager
// -----------------------------------------------------------------

VSStlsManager::VSStlsManager(const std::shared_ptr<const VSStlsInput> &in)
    : VSManager(),
      Stls(in, false),
      managerInPtr(in) {
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
      workers[idx] = std::make_unique<VSStlsWorker>(inTmp);
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
      grid_(in) {
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
