#include "vsstls_new.hpp"
#include "input.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;

// -----------------------------------------------------------------
// VSStls class
// -----------------------------------------------------------------

VSStlsNew::VSStlsNew(const VSStlsInput &in_)
  : VSBase(in_), Stls(in_.toStlsInput()), thermoProp(make_shared<ThermoProp>(in_)) {
  VSBase::thermoProp = thermoProp;
}

VSStlsNew::VSStlsNew(const VSStlsInput &in_, const ThermoProp &thermoProp_)
  : VSBase(in_), Stls(in_.toStlsInput(), false, false), thermoProp(make_shared<ThermoProp>(in_)) {
  VSBase::thermoProp = thermoProp;
  thermoProp->copyFreeEnergyIntegrand(thermoProp_);
}

double VSStlsNew::computeAlpha() {
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
  double numer = 2 * fxc - (1.0 / 6.0) * fxcrr + (4.0 / 3.0) * fxcr;
  double denom = uint + (1.0 / 3.0) * uintr;
  if (in.getDegeneracy() > 0.0) {
    numer += -(2.0 / 3.0) * fxctt - (2.0 / 3.0) * fxcrt + (1.0 / 3.0) * fxct;
    denom += (2.0 / 3.0) * uintt;
  }
  return numer / denom;
}

void VSStlsNew::updateSolution() {
  // Update the structural properties used for output
  slfc = thermoProp->getSlfc();
  ssf = thermoProp->getSsf();
}

void VSStlsNew::initScheme() {
  Rpa::init();
}

void VSStlsNew::initFreeEnergyIntegrand() {
  if (!thermoProp->isFreeEnergyIntegrandIncomplete()) { return; }
  if (verbose) {
    printf("Missing points in the free energy integrand: subcalls will be "
           "performed to collect the necessary data\n");
  }
  if (verbose) {
    printf("-----------------------------------------------------------------"
           "----------\n");
  }
  VSStlsInput inTmp = in;
  while (thermoProp->isFreeEnergyIntegrandIncomplete()) {
    const double rs = thermoProp->getFirstUnsolvedStatePoint();
    if (verbose) { printf("Subcall: solving VS scheme for rs = %.5f:\n", rs); }
    inTmp.setCoupling(rs);
    VSStlsNew scheme(inTmp, *thermoProp);
    scheme.compute();
    thermoProp->copyFreeEnergyIntegrand(*(scheme.thermoProp));
    if (verbose) {
      printf("Done\n");
      printf("-----------------------------------------------------------------"
             "----------\n");
    }
  }
}

// -----------------------------------------------------------------
// ThermoPropBase class
// -----------------------------------------------------------------

ThermoProp::ThermoProp(const VSStlsInput &in_): ThermoPropBase(in_), structProp(make_shared<StructPropNew>(in_)) {
  ThermoPropBase::structProp = structProp;
}


// -----------------------------------------------------------------
// StructProp class
// -----------------------------------------------------------------

StructPropNew::StructPropNew(const VSStlsInput &in_) : StructPropBase() {
  setupCSR(in_);
  setupCSRDependencies();
  for (const auto& c : csr) {
    StructPropBase::csr.push_back(c);
  }
}
						       
void StructPropNew::setupCSR(const VSStlsInput &in) {
  std::vector<VSStlsInput> inVector = setupCSRInput(in);
  for (const auto &inTmp : inVector) {
    csr.push_back(make_shared<StlsCSRNew>(inTmp));
  }
  assert(csr.size() == NPOINTS);
}

std::vector<VSStlsInput> StructPropNew::setupCSRInput(const VSStlsInput &in) {
  const double &drs = in.getCouplingResolution();
  const double &dTheta = in.getDegeneracyResolution();
  // If there is a risk of having negative state parameters, shift the
  // parameters so that rs - drs = 0 and/or theta - dtheta = 0
  const double rs = std::max(in.getCoupling(), drs);
  const double theta = std::max(in.getDegeneracy(), dTheta);
  // Setup objects
  std::vector<VSStlsInput> out;
  for (const double &thetaTmp : {theta - dTheta, theta, theta + dTheta}) {
    for (const double &rsTmp : {rs - drs, rs, rs + drs}) {
      VSStlsInput inTmp = in;
      inTmp.setDegeneracy(thetaTmp);
      inTmp.setCoupling(rsTmp);
      out.push_back(inTmp);
    }
  }
  return out;
}

void StructPropNew::doIterations() {
  const auto &in = csr[0]->getInput();
  const int maxIter = in.getNIter();
  const int ompThreads = in.getNThreads();
  const double minErr = in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  for (auto &c : csr) {
    c->initialGuess();
  }
  // Iteration to solve for the structural properties
  const bool useOMP = ompThreads > 1;
  while (counter < maxIter + 1 && err > minErr) {
// Compute new solution and error
#pragma omp parallel num_threads(ompThreads) if (useOMP)
    {
#pragma omp for
      for (auto &c : csr) {
        c->computeSsf();
        c->computeSlfcStls();
      }
#pragma omp for
      for (size_t i = 0; i < csr.size(); ++i) {
        auto &c = csr[i];
        c->computeSlfc();
        if (i == RS_THETA) { err = c->computeError(); }
        c->updateSolution();
      }
    }
    counter++;
  }
  if (verbose) {
    printf("Alpha = %.5e, Residual error "
           "(structural properties) = %.5e\n",
           csr[RS_THETA]->getAlpha(),
           err);
  }
}

// -----------------------------------------------------------------
// StlsCSR class
// -----------------------------------------------------------------

void StlsCSRNew::computeSlfcStls() {
  Stls::computeSlfc();
  *lfc = slfcNew;
}

void StlsCSRNew::computeSlfc() {
  // Check that alpha has been set to a value that is not the default
  assert(alpha != DEFAULT_ALPHA);
  // Derivative contributions
  const double &rs = in.getCoupling();
  // const double& theta = in.getDegeneracy();
  const double &theta = 0.0;
  const double &dx = in.getWaveVectorGridRes();
  const double &drs = in.getCouplingResolution();
  const double &dTheta = in.getDegeneracyResolution();
  const vector<double> &lfcData = *lfc;
  const vector<double> &rsUp = *lfcRs.up;
  const vector<double> &rsDown = *lfcRs.down;
  const vector<double> &thetaUp = *lfcTheta.up;
  const vector<double> &thetaDown = *lfcTheta.down;
  const double a_drs = alpha * rs / (6.0 * drs);
  const double a_dx = alpha / (6.0 * dx);
  const double a_dt = alpha * theta / (3.0 * dTheta);
  const size_t nx = wvg.size();
  // Wave-vector derivative
  slfcNew[0] -= a_dx * wvg[0] * getDerivative(lfc, 0, FORWARD);
  for (size_t i = 1; i < nx - 1; ++i) {
    slfcNew[i] -= a_dx * wvg[i] * getDerivative(lfc, i, CENTERED);
  }
  slfcNew[nx - 1] -= a_dx * wvg[nx - 1] * getDerivative(lfc, nx - 1, BACKWARD);
  // Coupling parameter contribution
  if (rs > 0.0) {
    for (size_t i = 0; i < nx; ++i) {
      slfcNew[i] -= a_drs * CSR::getDerivative(
                                lfcData[i], rsUp[i], rsDown[i], lfcRs.type);
    }
  }
  // Degeneracy parameter contribution
  if (theta > 0.0) {
    for (size_t i = 0; i < nx; ++i) {
      slfcNew[i] -=
          a_dt * CSR::getDerivative(
                     lfcData[i], thetaUp[i], thetaDown[i], lfcTheta.type);
    }
  }
}

double StlsCSRNew::getDerivative(const shared_ptr<vector<double>> &f,
				 const size_t &idx,
				 const Derivative &type) {
  const vector<double> fData = *f;
  switch (type) {
  case BACKWARD:
    assert(idx >= 2);
    return CSR::getDerivative(fData[idx], fData[idx - 1], fData[idx - 2], type);
    break;
  case CENTERED:
    assert(idx >= 1 && idx < fData.size() - 1);
    return CSR::getDerivative(fData[idx], fData[idx + 1], fData[idx - 1], type);
    break;
  case FORWARD:
    assert(idx < fData.size() - 2);
    return CSR::getDerivative(fData[idx], fData[idx + 1], fData[idx + 2], type);
    break;
  default:
    assert(false);
    return -1;
    break;
  }
}
