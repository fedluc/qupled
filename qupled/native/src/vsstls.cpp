#include "vsstls.hpp"
#include "input.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"
#include <fmt/core.h>

using namespace std;

// -----------------------------------------------------------------
// VSStls class
// -----------------------------------------------------------------

VSStls::VSStls(const VSStlsInput &in_)
    : VSBase(in_),
      Stls(in_, false, false),
      in(in_),
      thermoProp(make_shared<ThermoProp>(in_)) {
  VSBase::thermoProp = thermoProp;
}

VSStls::VSStls(const VSStlsInput &in_, const ThermoProp &thermoProp_)
    : VSBase(in_, false),
      Stls(in_, false, false),
      in(in_),
      thermoProp(make_shared<ThermoProp>(in_)) {
  VSBase::thermoProp = thermoProp;
  thermoProp->copyFreeEnergyIntegrand(thermoProp_);
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
  double numer = 2 * fxc - (1.0 / 6.0) * fxcrr + (4.0 / 3.0) * fxcr;
  double denom = uint + (1.0 / 3.0) * uintr;
  if (in.getDegeneracy() > 0.0) {
    numer += -(2.0 / 3.0) * fxctt - (2.0 / 3.0) * fxcrt + (1.0 / 3.0) * fxct;
    denom += (2.0 / 3.0) * uintt;
  }
  return numer / denom;
}

void VSStls::updateSolution() {
  // Update the structural properties used for output
  slfc = thermoProp->getSlfc();
  ssf = thermoProp->getSsf();
}

void VSStls::initScheme() { Rpa::init(); }

void VSStls::initFreeEnergyIntegrand() {
  if (!thermoProp->isFreeEnergyIntegrandIncomplete()) { return; }
  println("Missing points in the free energy integrand: subcalls will be "
          "performed to collect the necessary data");
  println("-----------------------------------------------------------------"
          "----------");
  VSStlsInput inTmp = in;
  while (thermoProp->isFreeEnergyIntegrandIncomplete()) {
    const double rs = thermoProp->getFirstUnsolvedStatePoint();
    println(fmt::format("Subcall: solving VS scheme for rs = {:.5f}", rs));
    inTmp.setCoupling(rs);
    VSStls scheme(inTmp, *thermoProp);
    scheme.compute();
    thermoProp->copyFreeEnergyIntegrand(*scheme.thermoProp);
    println("Done");
    println("-----------------------------------------------------------------"
            "----------");
  }
  println("Subcalls completed");
}

// -----------------------------------------------------------------
// ThermoPropBase class
// -----------------------------------------------------------------

ThermoProp::ThermoProp(const VSStlsInput &in_)
    : ThermoPropBase(in_, in_),
      structProp(make_shared<StructProp>(in_)) {
  ThermoPropBase::structProp = structProp;
}

// -----------------------------------------------------------------
// StructProp class
// -----------------------------------------------------------------

StructProp::StructProp(const VSStlsInput &in_)
    : Logger(MPIUtil::isRoot()),
      StructPropBase(),
      in(in_) {
  setupCSR();
  setupCSRDependencies();
}

void StructProp::setupCSR() {
  std::vector<VSStlsInput> inVector = setupCSRInput();
  for (const auto &inTmp : inVector) {
    csr.push_back(make_shared<StlsCSR>(inTmp));
  }
  for (const auto &c : csr) {
    StructPropBase::csr.push_back(c);
  }
}

std::vector<VSStlsInput> StructProp::setupCSRInput() {
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

void StructProp::doIterations() {
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
  println(fmt::format("Alpha = {:.5e}, Residual error "
                      "(structural properties) = {:.5e}",
                      csr[RS_THETA]->getAlpha(),
                      err));
}

// -----------------------------------------------------------------
// StlsCSR class
// -----------------------------------------------------------------

void StlsCSR::computeSlfcStls() {
  Stls::computeSlfc();
  *lfc = Vector2D(slfcNew);
}

void StlsCSR::computeSlfc() {
  Vector2D slfcDerivative = getDerivativeContribution();
  for (size_t i = 0; i < slfcNew.size(); ++i) {
    slfcNew[i] -= slfcDerivative(i);
  }
}
