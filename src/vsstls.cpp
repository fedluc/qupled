#include <omp.h>
#include "vsstls.hpp"

using namespace vecUtil;
using namespace stringUtil;
using namespace thermoUtil;
using namespace binUtil;

// -----------------------------------------------------------------
// StlsCSR class
// -----------------------------------------------------------------

void StlsCSR::setDerivativeData(std::vector<std::shared_ptr<StlsCSR>>& stlsVector,
				const size_t& thisIdx) {
  // Find the index that corresponds to this state point in the vector
  const double& rs = in.getCoupling();
  const double& theta = in.getDegeneracy();
  const double& otherRs = stlsVector[thisIdx]->in.getCoupling();
  const double& otherTheta = stlsVector[thisIdx]->in.getDegeneracy();  
  assert(otherRs == rs && otherTheta == theta);
  const int&  STENCIL = StructProp::STENCIL;
  // Set data for coupling parameter derivative
  if (thisIdx % STENCIL == 0) {
    // Forward difference for all state points with rs - drs
    setDrsData(stlsVector[thisIdx + 1]->slfc, stlsVector[thisIdx + 2]->slfc, FORWARD);
  }
  else if ( thisIdx % STENCIL == 2 ) {
    // Backward difference for all state points with rs + drs
    setDrsData(stlsVector[thisIdx - 1]->slfc, stlsVector[thisIdx - 2]->slfc, BACKWARD);
  }
  else {
    // Centered difference for all state points with rs
    setDrsData(stlsVector[thisIdx + 1]->slfc, stlsVector[thisIdx - 1]->slfc, CENTERED);
  }
  // Set data for degeneracy parameter derivative
  if (thisIdx/STENCIL == 0) {
    // Forward difference for all state points with theta - dtheta
    setDThetaData(stlsVector[thisIdx + 1]->slfc, stlsVector[thisIdx + 2]->slfc, FORWARD);
  }
  else if ( thisIdx/STENCIL == STENCIL - 1 ) {
    // Backward difference for all state points with theta - dtheta
    setDThetaData(stlsVector[thisIdx - 1]->slfc, stlsVector[thisIdx - 2]->slfc, BACKWARD);
  }
  else {
    // Centered difference for all state points with theta
    setDThetaData(stlsVector[thisIdx + 1]->slfc, stlsVector[thisIdx - 1]->slfc, CENTERED);
  }
}

void StlsCSR::setDrsData(vector<double> &slfcStlsRsUp,
			 vector<double> &slfcStlsRsDown,
			 const Derivative &dTypeRs) {
  this->slfcStlsRsUp = &slfcStlsRsUp;
  this->slfcStlsRsDown = &slfcStlsRsDown;
  this->dTypeRs = dTypeRs;
}

void StlsCSR::setDThetaData(vector<double> &slfcStlsThetaUp,
			    vector<double> &slfcStlsThetaDown,
			    const Derivative &dTypeTheta) {
  this->slfcStlsThetaUp = &slfcStlsThetaUp;
  this->slfcStlsThetaDown = &slfcStlsThetaDown;
  this->dTypeTheta = dTypeTheta;
}

double StlsCSR::getDerivative(const vector<double>& f,
			      const size_t& idx,
			      const Derivative& type) {
  switch(type) {
  case BACKWARD:
    assert(idx - 2 >= 0);
    return getDerivative(f[idx], f[idx - 1], f[idx - 2], type);
    break;
  case CENTERED:
    assert(idx - 1 >= 0 && idx + 1 < f.size());
    return getDerivative(f[idx], f[idx + 1], f[idx - 1], type);
    break;
  case FORWARD:
    assert(idx + 2 < f.size());
    return getDerivative(f[idx], f[idx + 1], f[idx + 2], type);
    break;
  }
}

double StlsCSR::getDerivative(const double& f0,
			      const double& f1,
			      const double& f2,
			      const Derivative& type) {
  switch(type) {
  case BACKWARD:
    return 3.0 * f0 - 4.0 * f1 + f2; break;
  case CENTERED:
    return f1 - f2; break;
  case FORWARD:
    return -getDerivative(f0, f1, f2, BACKWARD); break;
  }
}

void StlsCSR::computeSlfcStls() {
  Stls::computeSlfc();
  slfcStls = slfc;
}

void StlsCSR::computeSlfc() {
  // Derivative contributions
  const double& rs = in.getCoupling();
  const double& theta = in.getDegeneracy();
  const double& dx = in.getWaveVectorGridRes();
  const double& drs = in.getCouplingResolution();
  const double& dTheta = in.getDegeneracyResolution();
  const vector<double>& rsUp = *slfcStlsRsUp;
  const vector<double>& rsDown = *slfcStlsRsDown;
  const vector<double>& thetaUp = *slfcStlsThetaUp;
  const vector<double>& thetaDown = *slfcStlsThetaDown;
  const double a_drs = alpha * rs / (6.0 * drs);
  const double a_dx = alpha/(6.0 * dx);
  const double a_dtheta = alpha * theta / (3.0 * dTheta);
  for (size_t i = 0; i < wvg.size(); ++i) {
    Derivative dType = CENTERED;
    if (i == 0) { dType = FORWARD; }
    else if (i == wvg.size() - 1) { dType = BACKWARD; }
    // Wave-vector derivative
    slfc[i] -= a_dx * wvg[i] * getDerivative(slfcStls, i, dType);
    // Coupling parameter derivative
    if (rs > 0) {
      slfc[i] -= a_drs * getDerivative(slfcStls[i], rsUp[i], rsDown[i], dTypeRs);
    }
    // Degeneracy parameter derivative
    if (theta > 0) { 
      slfc[i] -= a_dtheta * getDerivative(slfcStls[i], thetaUp[i], thetaDown[i], dTypeTheta);
    }
  }
}

double StlsCSR::getFreeEnergyIntegrand() const {
  return computeInternalEnergy(wvg, ssf, 1.0);
}

// -----------------------------------------------------------------
// StructProp class
// -----------------------------------------------------------------

StructProp::StructProp(const VSStlsInput &in_) : in(in_), stlsIsInitialized(false) {
  VSStlsInput inTmp = in;
  const double& drs = inTmp.getCouplingResolution();
  const double& dTheta = inTmp.getDegeneracyResolution();
  // If there is a risk of having negative state parameters, shift the
  // parameters so that rs - drs = 0 and/or theta - dtheta = 0
  if (inTmp.getCoupling() <= drs) { inTmp.setCoupling(2.0 * drs); }
  if (inTmp.getDegeneracy() < dTheta) { inTmp.setDegeneracy(dTheta); }
  double rs = inTmp.getCoupling();
  double theta = inTmp.getDegeneracy();
  // Setup objects  
  for (int i = -1; i < 2; ++i){
    for (int j = -1; j < 2; ++j) {
      inTmp.setDegeneracy(theta + i * dTheta);
      inTmp.setCoupling(rs + j * drs);
      stls.push_back(std::make_shared<StlsCSR>(inTmp));
    }
  }
  assert(stls.size() == NPOINTS);
  // Setup derivative dependency in the StlsCSR objects
  for (size_t i = 0; i < stls.size(); ++i) {
    assert(stls[i]);
    stls[i]->setDerivativeData(stls, i);
  }
}

// Add a public call to compute where we do the initializations
int StructProp::compute() {
  try {
    if (!stlsIsInitialized) {
      for (auto& s : stls) { s->init(); }
      stlsIsInitialized = true;
    }
    doIterations();
    return 0;
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    return 1;
  }
}

void StructProp::doIterations() {
  const int maxIter = in.getNIter();
  const double minErr = in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  for (auto& s : stls) { s->initialGuess(); }
  // Iteration to solve for the structural properties
  while (counter < maxIter+1 && err > minErr ) {
    // Compute new solution and error
    for (auto& s : stls) { s->computeSsf(); }
    for (auto& s : stls) { s->computeSlfcStls(); }
    for (auto& s : stls) { s->computeSlfc(); }
    // Update error and diagnostic
    counter++;
    err = 0;
    for (auto& s : stls) { err +=  s->computeError(); }
    // Update solution
    for (auto& s : stls) { s->updateSolution(); }
  }
}

void StructProp::setAlpha(const double& alpha) {
  for (auto& s : stls) { s->setAlpha(alpha); }
}

vector<double> StructProp::getFreeEnergyIntegrand(const double& theta) {
  vector<double> out;
  for (auto& s : stls) {
    if (s->in.getDegeneracy() == theta) { out.push_back(s->getFreeEnergyIntegrand()); }
  }
  assert(out.size() == STENCIL);
  return out;
}

const StlsCSR& StructProp::getStls(const double& rs,
				   const double& theta) const {
  for (auto& s : stls) {
    if (s->in.getCoupling() == rs && s->in.getDegeneracy() == theta) {
      return *s;
    }
  }
  assert(false);
  return *stls[0];
}

// -----------------------------------------------------------------
// VSStls class
// -----------------------------------------------------------------

VSStls::VSStls(const VSStlsInput &in_) : StlsBase(in_), in(in_), structProp(in_) {
  const double& rs = in.getCoupling();
  const double& drs = in.getCouplingResolution();
  // Build integration grid
  rsGrid.push_back(0.0);
  const double rsMax = rs + drs;
  while(rsGrid.back() < rsMax){
    rsGrid.push_back(rsGrid.back() + drs);
  }
  // Resize the vector to store the free energy integrand;
  freeEnergyIntegrand.resize(StructProp::STENCIL);
  for (auto& f : freeEnergyIntegrand) {
    f.resize(rsGrid.size());
    fill(f, 0.0);
  }
}

int VSStls::compute() {
  try {
    init();  
    // if (verbose) cout << "Structural properties calculation ..." << endl;
    doIterations();
    // if (verbose) cout << "Done" << endl;
    return 0;
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    return 1;
  }
}

// Initialization
void VSStls::init() {
  computeFixedFreeEnergy();
}

// Compute fixed component of the free energy
void VSStls::computeFixedFreeEnergy() {
  // Compute the free energy from file passed in input, skipped for now
  return; 
}

// stls iterations
void VSStls::doIterations() {
  const int maxIter = in.getNIter();
  const int outIter = in.getOutIter();
  const double minErr = in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  initialGuess();
  while (counter < maxIter+1 && err > minErr ) {
    // Start timing
    double tic = omp_get_wtime();
    // Update static structure factor
    computeAlpha();
    // Update diagnostic
    counter++;
    err = computeError();
    // Update solution
    updateSolution();
    // Write output
    // if (counter % outIter == 0 && writeFiles) { writeRecovery(); }
    // End timing
    double toc = omp_get_wtime();
    // Print diagnostic
    printf("--- iteration %d ---\n", counter);
    printf("Elapsed time: %f seconds\n", toc - tic);
    printf("Residual error: %.5e\n", err);
    printf("alpha (CSR): %.5e\n", alpha);
    fflush(stdout);
  }
}

void VSStls::initialGuess() {
  alpha = in.getAlpha();
}

void VSStls::computeAlpha() {
  const double drs = in.getCouplingResolution();
  const double dt = in.getDegeneracyResolution();
  const double rs = in.getCoupling();
  const double theta = in.getDegeneracy();
  const double rsm = rs - drs;
  const double rsp = rs + drs;
  const double rs2 = rs * rs;      
  // // Check if a finite temperature is used
  // const bool isFiniteTemperature = (in.getDegeneracy() > 0.0);
  // assert(isFiniteTemperature);
  // // Compute the free energy integrand
  // computeFreeEnergyIntegrand();
  // // Compute internal energy
  // const size_t nrs = rsGrid.size();
  // const double urst = freeEnergyIntegrand[1][nrs - 2] / rs;
  // // if (isFiniteTemperature) {
  // const double urstm = freeEnergyIntegrand[0][nrs - 2] / rs;
  // const double urstp = freeEnergyIntegrand[2][nrs - 2] / rs;
  // // }
  // // Compute free energy
  // const double frsmt = computeFreeEnergy(rsGrid, freeEnergyIntegrand[1], rsm);
  // const double frst = computeFreeEnergy(rsGrid, freeEnergyIntegrand[1], rs);
  // const double frspt = computeFreeEnergy(rsGrid, freeEnergyIntegrand[1], rsp);
  // // if (finite_temperature) {
  // const double frsmtm = computeFreeEnergy(rsGrid, freeEnergyIntegrand[0], rsm);
  // const double frstm = computeFreeEnergy(rsGrid, freeEnergyIntegrand[0], rs);
  // const double frsptm = computeFreeEnergy(rsGrid, freeEnergyIntegrand[0], rsp);
  // const double frsmtp = computeFreeEnergy(rsGrid, freeEnergyIntegrand[2], rsm);
  // const double frstp = computeFreeEnergy(rsGrid, freeEnergyIntegrand[2], rs);
  // const double frsptp = computeFreeEnergy(rsGrid, freeEnergyIntegrand[2], rsp);
  // // }
  // // Internal energy derivatives
  // const double dudrs = (freeEnergyIntegrand[1][nrs - 1] - freeEnergyIntegrand[1][nrs - 3])/(2.0 * drs) - urst;
  // // if (finite_temperature)
  // const double dudt =  (urstp - urstm)/(2.0*dt);
  // // Free energy derivatives
  // const double dfdrs = (frspt - frsmt)/(2.0*drs);
  // const double d2fdrs2 = (frspt -2.0*frst + frsmt)/(drs*drs);
  // // if (finite_temperature) {
  // const double dfdt = (frstp - frstm)/(2.0*dt);
  // const double d2fdt2 = (frstp -2.0*frst + frstm)/(dt*dt);
  // const double d2fdrsdt = (frsptp - frsmtp - frsptm + frsmtm)/(4.0*drs*dt);
  // //  }
  // // Parameter for the compressibility sum rule
  // double numer = 2.0*frst - (1.0/6.0) * rs2 * d2fdrs2 + (4.0/3.0) * rs * dfdrs;
  // double denom = urst + dudrs/3.0 ; // (1.0/3.0) * rs * dudrs; 
  // numer += -(2.0/3.0) * theta * theta * d2fdt2
  //   -(2.0/3.0) * theta * rs * d2fdrsdt
  //   +(1.0/3.0) * theta * dfdt;
  // denom += (2.0/3.0) * theta *dudt;
  // alphaNew =  numer/denom;
  // std::cerr << alphaNew << std::endl;
  const bool isFiniteTemperature = (in.getDegeneracy() > 0.0);
  assert(!isFiniteTemperature);
  // Compute the free energy integrand
  computeFreeEnergyIntegrand();
  const size_t nrs = rsGrid.size();
  // Free energy
  const double fxc = computeFreeEnergy(rsGrid, freeEnergyIntegrand[0], rs);
  const double fxcUp = computeFreeEnergy(rsGrid, freeEnergyIntegrand[0], rsp);
  const double fxcDown = computeFreeEnergy(rsGrid, freeEnergyIntegrand[0], rsm);
  const double dfxcdrs = (fxcUp - fxcDown) / (2.0 * drs);
  const double d2fxcdrs = (fxcUp - 2.0 * fxc + fxcDown) / (drs * drs);
  // Internal energy
  const double uint = freeEnergyIntegrand[0][nrs - 2] / rs;
  const double drsudrs = (freeEnergyIntegrand[0][nrs - 1]
			  - freeEnergyIntegrand[0][nrs - 3]) / (2.0 * drs) - uint;
  // Alpha
  const double numer = 2 * fxc - (1.0/6.0) * rs2 * d2fxcdrs + (4.0/3.0) * rs * dfxcdrs;
  const double denom = uint + (1.0/3.0) * rs * drsudrs;
  alphaNew = numer/denom;
}

double VSStls::computeError() {
  const double diff = alpha - alphaNew;
  return sqrt(diff*diff);
}

void VSStls::updateSolution() {
  // Update the free parameter
  const double aMix = in.getMixingParameterAlpha();
  alpha *= (1 - aMix);
  alpha += alphaNew * aMix;
  structProp.setAlpha(alpha);
  // Update the structural properties used for output
  const auto& stls = structProp.getStls(in.getCoupling(), in.getDegeneracy());
  wvg = stls.getWvg();
  idr = stls.getIdr();
  slfc = stls.getSlfc();
  ssf = stls.getSsf();
  ssfHF = stls.getSsfHF();
}

void VSStls::computeFreeEnergyIntegrand() {
  const double drs = in.getCouplingResolution();
  const double dTheta = in.getDegeneracyResolution();
  const double Theta = in.getDegeneracy();
  VSStlsInput inTmp = in;
  for (size_t i = 0; i < rsGrid.size(); ++i) {
    const double& rs = rsGrid[i];
    if (rs < in.getCoupling()) {
      if (rs == 0) {
	std::cout << "computeFreeEnergyIntegrand skipping rs = " << rs << std::endl;
      }
      else if (rs == drs) {
	std::cout << "computeFreeEnergyIntegrand solving rs = " << rs << std::endl;
	inTmp.setCoupling(rs);
	VSStls vsstlsTmp(inTmp);
	vsstlsTmp.compute();
      }
      else {
	// std::cout << "computeFreeEnergyIntegrand nothing to do for rs = " << rs << std::endl;
      }
    }
  }
  if (in.getCoupling() == drs) {
    std::cerr << "structProp.compute() for rs = " << in.getCoupling() << std::endl;
    structProp.compute();
    // freeEnergyIntegrand[0] = structProp.getFreeEnergyIntegrand(Theta - dTheta);
    // freeEnergyIntegrand[1] = structProp.getFreeEnergyIntegrand(Theta);
    // freeEnergyIntegrand[2] = structProp.getFreeEnergyIntegrand(Theta + dTheta);
    freeEnergyIntegrand[0] = structProp.getFreeEnergyIntegrand(Theta);
  }
}
