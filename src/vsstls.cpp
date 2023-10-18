#include <omp.h>
#include "vsstls.hpp"

using namespace vecUtil;
using namespace stringUtil;
using namespace thermoUtil;
using namespace binUtil;

// -----------------------------------------------------------------
// StlsCSR class
// -----------------------------------------------------------------

void StlsCSR::setDerivativeData(std::vector<std::unique_ptr<StlsCSR>>& stlsVector,
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
  else if ( (thisIdx - STENCIL + 1) % STENCIL == 0 ) {
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

void StlsCSR::setDrsData(vector<double> &slfcRsUp,
			 vector<double> &slfcRsDown,
			 const Derivative &dTypeRs) {
  this->slfcRsUp = &slfcRsUp;
  this->slfcRsDown = &slfcRsDown;
  this->dTypeRs = dTypeRs;
}

void StlsCSR::setDThetaData(vector<double> &slfcThetaUp,
			    vector<double> &slfcThetaDown,
			    const Derivative &dTypeTheta) {
  this->slfcThetaUp = &slfcThetaUp;
  this->slfcThetaDown = &slfcThetaDown;
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

double StlsCSR::doAction(const Action& action) {
  switch (action) {
  case INITIALIZE : Stls::init(); return -1;
  case GUESS: Stls::initialGuess(); return -1;
  case SOLUTION: Stls::computeSsf(); computeSlfc(); return -1;
  case ERROR: return Stls::computeError();
  case UPDATE: Stls::updateSolution(); return -1;
  default: return -1;
  }
}

void StlsCSR::computeSlfc() {
  // Stls contribution
  Stls::computeSlfcStls();
  // Derivative contributions
  const double& rs = in.getCoupling();
  const double& theta = in.getDegeneracy();
  const double& dx = in.getWaveVectorGridRes();
  const double& drs = in.getCouplingResolution();
  const double& dTheta = in.getDegeneracyResolution();
  const vector<double>& rsUp = *slfcRsUp;
  const vector<double>& rsDown = *slfcRsDown; 
  const vector<double>& thetaUp = *slfcThetaUp;
  const vector<double>& thetaDown = *slfcThetaDown;
  const double a_drs = alpha * rs / (6.0 * drs);
  const double a_dx = alpha/(6.0 * dx);
  const double a_dtheta = alpha * theta / (3.0 * dTheta);
  for (size_t i = 0; i < wvg.size(); ++i) {
    Derivative dType = CENTERED;
    if (i == 0) { dType = FORWARD; }
    else if (i == wvg.size() - 1) { dType = BACKWARD; }
    // Wave-vector derivative
    slfc[i] -= a_dx * wvg[i] * getDerivative(slfc, i, dType);
    // Coupling parameter derivative
    slfc[i] -= a_drs * getDerivative(slfc[i], rsUp[i], rsDown[i], dTypeRs);
    // Degeneracy parameter derivative
    slfc[i] -= a_dtheta * getDerivative(slfc[i], thetaUp[i], thetaDown[i], dTypeTheta);
  }
}


// -----------------------------------------------------------------
// StructProp class
// -----------------------------------------------------------------

StructProp::StructProp(const VSStlsInput &in_) : in(in_), stls(NPOINTS) {
  const double& drs = in.getCouplingResolution();
  const double& dTheta = in.getDegeneracyResolution();
  double rs = in.getCoupling();
  double theta = in.getDegeneracy();
  // If there is a risk of having negative state parameters, shift the
  // parameters so that rs - drs = 0 and/or theta - dtheta = 0
  if (rs < drs) { rs = drs; }
  if (theta < dTheta) { theta = dTheta; }
  VSStlsInput inTmp = in;
  size_t cnt = 0;
  // Setup objects  
  for (int i = -1; i < 1; ++i){
    for (int j = -1; j < 1; ++j) {
      inTmp.setDegeneracy(rs + j * drs);
      inTmp.setCoupling(theta + i * dTheta);
      assert(cnt < NPOINTS);
      stls[cnt] = std::make_unique<StlsCSR>(inTmp);
      cnt++;
    }
  }
  // Setup derivative dependency in the StlsCSR objects
  for (size_t i = 0; i < stls.size(); ++i) {
    stls[i]->setDerivativeData(stls, i);
  }
}

// Add a public call to compute where we do the initializations
int StructProp::compute() {
  try {
    for (auto& s : stls) { s->doAction(StlsCSR::Action::INITIALIZE); }  
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
  for (auto& s : stls) { s->doAction(StlsCSR::Action::GUESS); }
  // Iteration to solve for the structural properties
  while (counter < maxIter+1 && err > minErr ) {
    // Compute new solution and error
    for (auto& s : stls) { s->doAction(StlsCSR::Action::SOLUTION); }
    // Update diagnostic
    counter++;
    for (auto& s : stls) { s->doAction(StlsCSR::Action::ERROR); }
    // Update solution
    for (auto& s : stls) { s->doAction(StlsCSR::Action::UPDATE); }
  }
}

void StructProp::setAlpha(const double& alpha) {
  for (auto& s : stls) { s->setAlpha(alpha); }
}

// -----------------------------------------------------------------
// ThermoProp class
// -----------------------------------------------------------------

ThermoProp::ThermoProp(const VSStlsInput &in) : prop(NPOINTS) {
  const double& rs = in.getCoupling();
  const double& drs = in.getCouplingResolution();
  // Build integration grid
  grid.push_back(0.0);
  const double rsMax = rs + 2.0*drs; // 2.0*drs is added to circumvent GSL errors
  while(grid.back() < rsMax){
    grid.push_back(grid.back() + drs);
  }
  // Resize vectors used to store the result of the integration
  const size_t sz = grid.size();
  for (auto& p : prop) {
    p.resize(sz);
  }
}

// -----------------------------------------------------------------
// VSStls class
// -----------------------------------------------------------------

int VSStls::compute(){
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
void VSStls::init(){
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
  // Check if a finite temperature is used
  const bool isFiniteTemperature = (in.getDegeneracy() > 0.0);
  // Compute the free energy integrand
  computeFreeEnergyIntegrand();
}

double VSStls::computeError() {
  const double diff = alpha - alphaNew;
  return sqrt(diff*diff);
}

void VSStls::updateSolution() {
  const double aMix = in.getMixingParameterAlpha();
  alpha *= (1 - aMix);
  alpha += alphaNew * aMix;
  structProp.setAlpha(alpha);
}

void computeFreeEnergyIntegrand() {
  std::cerr << "computeFreeEnergyIntegrand in not implemented" << std::endl;
}
