#include "vsstls.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;
using namespace MPIUtil;

// -----------------------------------------------------------------
// VSStls class
// -----------------------------------------------------------------

VSStls::VSStls(const std::shared_ptr<const VSStlsInput> &in_)
    : VSBase(),
      Stls(in_, false),
      thermoProp(make_shared<ThermoProp>(in_)) {
  if (in_->getDimension() == dimensionsUtil::Dimension::D2) {
    throwError("2D calculations are not implemented for this scheme.");
  }
  VSBase::thermoProp = thermoProp;
}

double VSStls::computeAlpha() {
  // Compute the free energy integrand
  thermoProp->compute();
  // Free energy
  const vector<double> freeEnergyData = thermoProp->getFreeEnergyData();
  const double &fxc = freeEnergyData[0];
  const double &fxcr = freeEnergyData[1];
  const double &fxcrr = freeEnergyData[2];
  const double &fxct = freeEnergyData[3];
  const double &fxctt = freeEnergyData[4];
  const double &fxcrt = freeEnergyData[5];
  // Internal energy
  const vector<double> internalEnergyData = thermoProp->getInternalEnergyData();
  const double &uint = internalEnergyData[0];
  const double &uintr = internalEnergyData[1];
  const double &uintt = internalEnergyData[2];
  // Alpha
  double numer = 2.0 * fxc + (4.0 / 3.0) * fxcr - (1.0 / 6.0) * fxcrr
                 - (2.0 / 3.0) * (fxctt + fxcrt) + (1.0 / 3.0) * fxct;
  double denom = uint + (1.0 / 3.0) * uintr + (2.0 / 3.0) * uintt;
  return numer / denom;
}

void VSStls::updateSolution() {
  // Update the structural properties used for output
  lfc = thermoProp->getLfc();
  ssf = thermoProp->getSsf();
}

void VSStls::init() { Rpa::init(); }

// -----------------------------------------------------------------
// ThermoProp class
// -----------------------------------------------------------------

ThermoProp::ThermoProp(const std::shared_ptr<const VSStlsInput> &in_)
    : ThermoPropBase(in_),
      structProp(make_shared<StructProp>(in_)) {
  ThermoPropBase::structProp = structProp;
}

// -----------------------------------------------------------------
// StructProp class
// -----------------------------------------------------------------

StructProp::StructProp(const std::shared_ptr<const VSStlsInput> &in_)
    : StructPropBase(in_),
      csr(make_shared<StlsCSRNew>(in_, true)) {
  StructPropBase::csr = csr;
}

// -----------------------------------------------------------------
// StlsCSR class
// -----------------------------------------------------------------

StlsCSRNew::StlsCSRNew(const std::shared_ptr<const VSStlsInput> &in_,
                       const bool isMaster_)
    : CSRNew(isMaster_),
      Stls(in_, false) {
  if (isManager) { setupAuxiliaryStatePoints(*in_); }
}

int StlsCSRNew::compute() { return Stls::compute(); }

void StlsCSRNew::setupAuxiliaryStatePoints(const VSStlsInput &in) {
  const double &drs = in.getCouplingResolution();
  const double &dTheta = in.getDegeneracyResolution();
  // If there is a risk of having negative state parameters, shift the
  // parameters so that rs - drs = 0 and/or theta - dtheta = 0
  const double rs = std::max(in.getCoupling(), drs);
  const double theta = std::max(in.getDegeneracy(), dTheta);
  // Setup auxiliary state points
  for (const double &thetaTmp : {theta - dTheta, theta, theta + dTheta}) {
    for (const double &rsTmp : {rs - drs, rs, rs + drs}) {
      std::shared_ptr<VSStlsInput> inTmp = std::make_shared<VSStlsInput>(in);
      inTmp->setDegeneracy(thetaTmp);
      inTmp->setCoupling(rsTmp);
      workers.push_back(make_shared<StlsCSRNew>(inTmp, false));
    }
  }
  assert(workers.size() == NRS * NTHETA);
  setupDerivativeData();
}
