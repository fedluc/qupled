#include <omp.h>
#include "vsstls.hpp"

using namespace numUtil;
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
    setDrsData(stlsVector[thisIdx + 1]->slfcStls, stlsVector[thisIdx + 2]->slfcStls, FORWARD);
  }
  else if ( thisIdx % STENCIL == 2 ) {
    // Backward difference for all state points with rs + drs
    setDrsData(stlsVector[thisIdx - 1]->slfcStls, stlsVector[thisIdx - 2]->slfcStls, BACKWARD);
  }
  else {
    // Centered difference for all state points with rs
    setDrsData(stlsVector[thisIdx + 1]->slfcStls, stlsVector[thisIdx - 1]->slfcStls, CENTERED);
  }
  // Set data for degeneracy parameter derivative
  if (thisIdx/STENCIL == 0) {
    // Forward difference for all state points with theta - dtheta
    setDThetaData(stlsVector[thisIdx + 1]->slfcStls, stlsVector[thisIdx + 2]->slfcStls, FORWARD);
  }
  else if ( thisIdx/STENCIL == STENCIL - 1 ) {
    // Backward difference for all state points with theta - dtheta
    setDThetaData(stlsVector[thisIdx - 1]->slfcStls, stlsVector[thisIdx - 2]->slfcStls, BACKWARD);
  }
  else {
    // Centered difference for all state points with theta
    setDThetaData(stlsVector[thisIdx + 1]->slfcStls, stlsVector[thisIdx - 1]->slfcStls, CENTERED);
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
  //const double& theta = in.getDegeneracy();
  const double& theta = 0.0;
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
  if (inTmp.getCoupling() < drs) { inTmp.setCoupling(drs); }
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

VSStls::VSStls(const VSStlsInput &in_) : StlsBase(in_), in(in_),
					 structProp(in_), verbose(true) {
  const double& rs = in.getCoupling();
  const double& drs = in.getCouplingResolution();
  // Build integration grid
  rsGrid = std::make_shared<std::vector<double>>();
  rsGrid->push_back(0.0);
  const double rsMax = rs + drs;
  while(!equalTol(rsGrid->back(), rsMax)){
    rsGrid->push_back(rsGrid->back() + drs);
  }
  // Resize the vector to store the free energy integrand;
  fxcIntegrand = std::make_shared<doubleVector>();
  fxcIntegrand->resize(StructProp::STENCIL);
  const size_t nrs = rsGrid->size();
  for (auto& f : *fxcIntegrand) {
    f.resize(nrs);
    fill(f, Inf);
  }
}

VSStls::VSStls(const VSStlsInput &in_,
	       std::shared_ptr<std::vector<double>> &rsGrid_,
	       std::shared_ptr<doubleVector> &fxcIntegrand_)
  : StlsBase(in_), in(in_), structProp(in_), rsGrid(rsGrid_),
    fxcIntegrand(fxcIntegrand_), verbose(false) { ; }

int VSStls::compute() {
  try {
    init();
    if (verbose) cout << "Free parameter calculation ..." << endl;
    doIterations();
    if (verbose) cout << "Done" << endl;
    return 0;
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    return 1;
  }
}

void VSStls::init() {
  const auto& fxcIntegrandIn = in.getFreeEnergyIntegrand();
  const size_t nrs = rsGrid->size();
  const size_t nrsIn = fxcIntegrandIn.grid.size();
  if (nrsIn == 0) {
    return;
  }
  const Interpolator1D fxci(fxcIntegrandIn.grid, fxcIntegrandIn.integrand);
  const double rsMaxi = fxcIntegrandIn.grid.back();
  for (size_t i = 0; i < nrs; ++i) {
    const double& rs = rsGrid->at(i);
    if (rs > rsMaxi) { return; }
    fxcIntegrand->at(0)[i] = fxci.eval(rs);
  }
}

// stls iterations
void VSStls::doIterations() {
  const int maxIter = in.getNIter();
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
    // End timing
    double toc = omp_get_wtime();
    // Print diagnostic
    if (verbose) {
      printf("--- iteration %d ---\n", counter);
      printf("Elapsed time: %f seconds\n", toc - tic);
      printf("Residual error: %.5e\n", err);
      printf("alpha (CSR): %.5e\n", alpha);
      fflush(stdout);
    }
  }
}

void VSStls::initialGuess() {
  alpha = in.getAlpha();
}

void VSStls::computeAlpha() {
  const double drs = in.getCouplingResolution();
  const double rs = in.getCoupling();
  const double rsm = rs - drs;
  const double rsp = rs + drs;
  // Compute the free energy integrand
  computeFreeEnergyIntegrand();
  // Free energy
  const double fxc = computeFreeEnergy(rs, true);
  const double fxcThis = computeFreeEnergy(rs, false);
  const double fxcUp = computeFreeEnergy(rsp, false); // rs2 * fxc
  const double fxcDown = computeFreeEnergy(rsm, false); // rs2 * fxc
  const double dfxc_drs = (fxcUp - fxcDown) / (2.0 * drs * rs) - 2.0 * fxc;
  const double d2fxc_drs = (fxcUp - 2.0 * fxcThis + fxcDown) / (drs * drs) - 2.0 * fxc - 4.0 * dfxc_drs;
  // Internal energy
  const auto& fxci = fxcIntegrand->at(0);
  const size_t nrs = rsGridLocal.size();
  const double uint = fxci[nrs - 2] / rs;
  const double uUp = fxci[nrs - 1]; // rs * uint
  const double uDown = fxci[nrs - 3]; // rs * uint
  const double du_drs = (uUp - uDown) / (2.0 * drs) - uint;
  // Alpha
  const double numer = 2 * fxc - (1.0/6.0) * d2fxc_drs + (4.0/3.0) * dfxc_drs;
  const double denom =  uint + (1.0/3.0) * du_drs;
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
  if (!structProp.isComputed()) { structProp.compute(); }
  const auto& stls = structProp.getStls(in.getCoupling(), in.getDegeneracy());
  wvg = stls.getWvg();
  idr = stls.getIdr();
  slfc = stls.getSlfc();
  ssf = stls.getSsf();
  ssfHF = stls.getSsfHF();
}

void VSStls::computeFreeEnergyIntegrand() {
  // Recursive calls to solve the VS-STLS scheme for all state points
  // with coupling parameter smaller than rs
  const double Theta = in.getDegeneracy();
  const double nrs = rsGrid->size();
  VSStlsInput inTmp = in;
  for (size_t i = 0; i < nrs; ++i) {
    const double rs = rsGrid->at(i);
    if (equalTol(rs, in.getCoupling())) {
      structProp.compute();
      fxcIntegrand->at(0)[i-1] = structProp.getFreeEnergyIntegrand(Theta)[0];
      fxcIntegrand->at(0)[i] = structProp.getFreeEnergyIntegrand(Theta)[1];
      fxcIntegrand->at(0)[i+1] = structProp.getFreeEnergyIntegrand(Theta)[2];
      i += (i < nrs - 4) ? 3 : 0;
    }
    else if (rs < in.getCoupling()) {
      if (rs == 0 || fxcIntegrand->at(0)[i] != Inf) {
	continue;
      }
      printf("Free energy integrand calculation, solving VS-STLS scheme for rs = %.5f: ", rs);
      inTmp.setCoupling(rs);
      VSStls vsstlsTmp(inTmp, rsGrid, fxcIntegrand);
      vsstlsTmp.compute();
      if (verbose) printf("Done\n");
    }
    else {
      break;
    }
  }
}

double VSStls::computeFreeEnergy(const double& rs,
				 const bool& normalize) {
  // Construct a local version of the free energy integrand that
  // does not contain infinities. i.e. that it contains only
  // the values that have been computed
  const size_t nrs = rsGrid->size();
  const auto& fxci = fxcIntegrand->at(0);
  const auto& fxciBegin = fxci.begin();
  const auto& fxciEnd = fxci.end();
  if (rsGridLocal.size() == 0) {
    auto isInf = [&](const double& num) { return num == Inf; };
    const auto it = std::find_if(fxciBegin, fxciEnd, isInf);
    size_t nCopy = std::distance(fxciBegin, it);
    std::copy(rsGrid->begin(), rsGrid->begin() + nCopy, std::back_inserter(rsGridLocal));
    fxcIntegrandLocal.resize(rsGridLocal.size());
  }
  if (rsGridLocal.size() < nrs) {
    std::copy(fxciBegin, fxciBegin + rsGridLocal.size(), fxcIntegrandLocal.begin());
    return thermoUtil::computeFreeEnergy(rsGridLocal, fxcIntegrandLocal, rs, normalize);
  }
  return thermoUtil::computeFreeEnergy(*rsGrid, fxci, rs, normalize);
}
