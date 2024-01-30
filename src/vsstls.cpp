#include <omp.h>
#include "vsstls.hpp"

using namespace numUtil;
using namespace vecUtil;
using namespace thermoUtil;

// -----------------------------------------------------------------
// VSStls class
// -----------------------------------------------------------------

int VSStls::compute() {
  try {
    omp_set_num_threads(in.getNThreads());
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

// stls iterations
void VSStls::doIterations() {
  auto func = [this](double alphaTmp)->double{return alphaDifference(alphaTmp);};
  SecantSolver rsol(in.getErrMinAlpha(), in.getNIterAlpha());
  rsol.solve(func, in.getAlphaGuess());
  if (!rsol.success()) {
    throw runtime_error("VSStls: the root solver did not converge to the desired accuracy.");
  }
  alpha = rsol.getSolution();
  if (verbose) { cout << "Free parameter = " << alpha << endl; }
  updateSolution();
}

double VSStls::alphaDifference(const double& alphaTmp) {
  alpha = alphaTmp;
  thermoProp.setAlpha(alpha);
  const double alphaTheoretical = computeAlpha();
  return alpha - alphaTheoretical;
}

double VSStls::computeAlpha() {
  // Compute the free energy integrand
  thermoProp.compute(in);
  // Free energy
  const vector<double> freeEnergyData = thermoProp.getFreeEnergyData();
  const double& fxc = freeEnergyData[0];
  const double& fxcr = freeEnergyData[1];
  const double& fxcrr = freeEnergyData[2];
  const double& fxct = freeEnergyData[3];
  const double& fxctt = freeEnergyData[4];
  const double& fxcrt = freeEnergyData[5];
  // Internal energy
  const vector<double> internalEnergyData = thermoProp.getInternalEnergyData();
  const double& uint = internalEnergyData[0];
  const double& uintr = internalEnergyData[1];
  const double& uintt = internalEnergyData[2];
  // Alpha
  double numer = 2 * fxc - (1.0/6.0) * fxcrr + (4.0/3.0) * fxcr;
  double denom =  uint + (1.0/3.0) * uintr;
  if (in.getDegeneracy() > 0.0) {
    numer += - (2.0/3.0) * fxctt
             - (2.0/3.0) * fxcrt
             + (1.0/3.0) * fxct;
    denom += (2.0/3.0) * uintt;
  }
  return numer/denom;
}

void VSStls::updateSolution() {
  // Update the structural properties used for output
  const auto& stls = thermoProp.getStructProp();
  slfc = stls.getSlfc();
  ssf = stls.getSsf();
}


vector<vector<double>> VSStls::getFreeEnergyIntegrand() const {
  return thermoProp.getFreeEnergyIntegrand();
}

vector<double> VSStls::getFreeEnergyGrid() const {
  return thermoProp.getFreeEnergyGrid();
}

// -----------------------------------------------------------------
// ThermoProp class
// -----------------------------------------------------------------

ThermoProp::ThermoProp(const VSStlsInput &in) : structProp(in) {
  const double& rs = in.getCoupling();
  const double& drs = in.getCouplingResolution();
  // Check if we are solving for particular state points
  isZeroCoupling = (rs == 0.0);
  isZeroDegeneracy = (in.getDegeneracy() == 0.0);
  // Build integration grid
  rsGrid.push_back(0.0);
  const double rsMax = rs + drs;
  while(!equalTol(rsGrid.back(), rsMax)){
    rsGrid. push_back(rsGrid.back() + drs);
  }
  // Initialize the free energy integrand
  fxcIntegrand.resize(NPOINTS);
  const size_t nrs = rsGrid.size();
  for (auto& f : fxcIntegrand) {
    f.resize(nrs);
    fill(f, Inf);
  }
  // Fill the free energy integrand if passed in input
  const auto& fxciData = in.getFreeEnergyIntegrand();
  if (!fxciData.grid.empty()) {
    for (const auto& theta : {Idx::THETA_DOWN, Idx::THETA, Idx::THETA_UP}) {
      const Interpolator1D itp(fxciData.grid, fxciData.integrand[theta]);
      const double rsMaxi = fxciData.grid.back();
      for (size_t i = 0; i < nrs; ++i) {
	const double& rs = rsGrid[i];
	if (rs <= rsMaxi) { fxcIntegrand[theta][i] = itp.eval(rs); }
      }
    }
  }
}

ThermoProp::ThermoProp(const VSStlsInput &in,
		       const ThermoProp &other) : ThermoProp(in) {
  assert(other.rsGrid[1] - other.rsGrid[0] == rsGrid[1] - rsGrid[0]);
  const size_t nrs = rsGrid.size();
  const double rsMax = other.rsGrid.back();
  for (const auto& theta : {Idx::THETA_DOWN, Idx::THETA, Idx::THETA_UP}) {
    const auto& fxciBegin = fxcIntegrand[theta].begin();
    const auto& fxciEnd = fxcIntegrand[theta].end();
    const auto& it = std::find(fxciBegin, fxciEnd, Inf);
    size_t i =  std::distance(fxciBegin, it);
    while (i < nrs && rsGrid[i] < rsMax) {
      fxcIntegrand[theta][i] = other.fxcIntegrand[theta][i];
      ++i;
    }
  }
}

void ThermoProp::setAlpha(const double& alpha) {
  structProp.setAlpha(alpha);
}

void ThermoProp::compute(const VSStlsInput& in) {
  // Recursive calls to solve the VS-STLS scheme for all state points
  // with coupling parameter smaller than rs
  const double nrs = rsGrid.size();
  VSStlsInput inTmp = in;
  vector<double> fxciTmp(StructProp::NPOINTS);
  for (size_t i = 0; i < nrs; ++i) {
    const double& rs = rsGrid[i];
    if (equalTol(rs, in.getCoupling())) {
      structProp.compute();
      fxciTmp = structProp.getFreeEnergyIntegrand();
    }
    else if (rs < in.getCoupling()) {
      if (rs == 0.0 || fxcIntegrand[THETA][i] != Inf) { continue; }
      printf("Free energy integrand calculation, solving VS-STLS scheme for rs = %.5f:\n", rs);
      inTmp.setCoupling(rs);
      VSStls vsstlsTmp(inTmp, *this);
      vsstlsTmp.compute();
      fxciTmp = vsstlsTmp.getThermoProp().structProp.getFreeEnergyIntegrand();
      printf("Done\n");
      printf("---------------------------------------------------------------------------\n");
    }
    else {
      break;
    }
    fxcIntegrand[THETA_DOWN][i-1] = fxciTmp[SIdx::RS_DOWN_THETA_DOWN];
    fxcIntegrand[THETA_DOWN][i]   = fxciTmp[SIdx::RS_THETA_DOWN];
    fxcIntegrand[THETA_DOWN][i+1] = fxciTmp[SIdx::RS_UP_THETA_DOWN];
    fxcIntegrand[THETA][i-1]      = fxciTmp[SIdx::RS_DOWN_THETA];
    fxcIntegrand[THETA][i]        = fxciTmp[SIdx::RS_THETA];
    fxcIntegrand[THETA][i+1]      = fxciTmp[SIdx::RS_UP_THETA];
    fxcIntegrand[THETA_UP][i-1]   = fxciTmp[SIdx::RS_DOWN_THETA_UP];
    fxcIntegrand[THETA_UP][i]     = fxciTmp[SIdx::RS_THETA_UP];
    fxcIntegrand[THETA_UP][i+1]   = fxciTmp[SIdx::RS_UP_THETA_UP];
  }
}

const StlsCSR& ThermoProp::getStructProp() {
  if (!structProp.isComputed()) { structProp.compute(); }
  if (isZeroCoupling && isZeroDegeneracy) {
     return structProp.getStls(SIdx::RS_DOWN_THETA_DOWN); 
  }
  if (!isZeroCoupling && isZeroDegeneracy) {
     return structProp.getStls(SIdx::RS_THETA_DOWN); 
  }
  if (isZeroCoupling && !isZeroDegeneracy) {
     return structProp.getStls(SIdx::RS_DOWN_THETA); 
  }
  return structProp.getStls(SIdx::RS_THETA); 
}

double ThermoProp::computeFreeEnergy(const SIdx iStruct,
				     const bool normalize) const {
  Idx iThermo;
  switch (iStruct) {
  case SIdx::RS_DOWN_THETA_DOWN: case SIdx::RS_THETA_DOWN: case SIdx::RS_UP_THETA_DOWN:
    iThermo = THETA_DOWN; break;
  case SIdx::RS_DOWN_THETA: case SIdx::RS_THETA: case SIdx::RS_UP_THETA:
    iThermo = THETA; break;
  case SIdx::RS_DOWN_THETA_UP: case SIdx::RS_THETA_UP: case SIdx::RS_UP_THETA_UP:
    iThermo = THETA_UP; break;
  }
  const vector<double>& rs = structProp.getCouplingParameters();
  return thermoUtil::computeFreeEnergy(rsGrid, fxcIntegrand[iThermo], rs[iStruct], normalize);
}

vector<double> ThermoProp::getFreeEnergyData() const {
  const vector<double> rsVec = structProp.getCouplingParameters();
  const vector<double> thetaVec = structProp.getDegeneracyParameters();
  // Free energy
  const double fxc = computeFreeEnergy(SIdx::RS_THETA, true);
  // Free energy derivatives with respect to the coupling parameter
  double fxcr;
  double fxcrr;
  {
    const double rs = rsVec[SIdx::RS_THETA];
    const double drs = rsVec[SIdx::RS_UP_THETA] - rsVec[SIdx::RS_THETA];
    const double f0 = computeFreeEnergy(SIdx::RS_UP_THETA, false);
    const double f1 = computeFreeEnergy(SIdx::RS_THETA, false);
    const double f2 = computeFreeEnergy(SIdx::RS_DOWN_THETA, false);
    fxcr = (f0 - f2) / (2.0 * drs * rs) - 2.0 * fxc;
    fxcrr = (f0 - 2.0 * f1 + f2) / (drs * drs) - 2.0 * fxc - 4.0 * fxcr;
  }
  // Free energy derivatives with respect to the degeneracy parameter
  double fxct;
  double fxctt;
  {
    const double theta = thetaVec[SIdx::RS_THETA];
    const double theta2 = theta * theta;
    const double dt = thetaVec[SIdx::RS_THETA_UP] - thetaVec[SIdx::RS_THETA];
    const double f0 = computeFreeEnergy(SIdx::RS_THETA_UP, true);
    const double f1 = computeFreeEnergy(SIdx::RS_THETA_DOWN, true);
    fxct = theta * (f0 - f1) / (2.0 * dt);
    fxctt = theta2 * (f0 - 2.0 * fxc + f1) / (dt * dt);
  }
  // Free energy mixed derivatives
  double fxcrt;
  {
    const double t_rs = thetaVec[SIdx::RS_THETA] / rsVec[SIdx::RS_THETA];
    const double drs = rsVec[SIdx::RS_UP_THETA] - rsVec[SIdx::RS_THETA];
    const double dt = thetaVec[SIdx::RS_THETA_UP] - thetaVec[SIdx::RS_THETA];
    const double f0 = computeFreeEnergy(SIdx::RS_UP_THETA_UP, false);
    const double f1 = computeFreeEnergy(SIdx::RS_UP_THETA_DOWN, false);
    const double f2 = computeFreeEnergy(SIdx::RS_DOWN_THETA_UP, false);
    const double f3 = computeFreeEnergy(SIdx::RS_DOWN_THETA_DOWN, false);
    fxcrt = t_rs * (f0 - f1 - f2 + f3) / (4.0 * drs * dt) - 2.0 * fxct;
  }
  return vector<double>({fxc, fxcr, fxcrr, fxct, fxctt, fxcrt});
}

vector<double> ThermoProp::getInternalEnergyData() const {
  // Internal energy
  const double u = structProp.getInternalEnergy()[SIdx::RS_THETA];
  // Internal energy derivative with respect to the coupling parameter
  double ur;
  {
    const vector<double> rs = structProp.getCouplingParameters(); 
    const double drs = rs[SIdx::RS_UP_THETA] - rs[SIdx::RS_THETA];
    const vector<double> rsu = structProp.getFreeEnergyIntegrand();
    const double& u0 = rsu[SIdx::RS_UP_THETA];
    const double& u1 = rsu[SIdx::RS_DOWN_THETA];
    ur = (u0 - u1) / (2.0 * drs) - u;
  }
  // Internal energy derivative with respect to the degeneracy parameter
  double ut;
  {
    const vector<double> theta = structProp.getDegeneracyParameters();
    const double dt = theta[SIdx::RS_THETA_UP] - theta[SIdx::RS_THETA];
    const double u0 = structProp.getInternalEnergy()[SIdx::RS_THETA_UP];
    const double u1 = structProp.getInternalEnergy()[SIdx::RS_THETA_DOWN];
    ut = theta[SIdx::RS_THETA] * (u0 - u1) / (2.0 * dt);
  }
  return vector<double>({u, ur, ut});
}


// -----------------------------------------------------------------
// StructProp class
// -----------------------------------------------------------------

StructProp::StructProp(const VSStlsInput &in) : computed(false), outVector(NPOINTS) {
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
  for (const double& thetaTmp : {theta - dTheta, theta, theta + dTheta}) {
    for (const double& rsTmp : {rs - drs, rs, rs + drs}){
      inTmp.setDegeneracy(thetaTmp);
      inTmp.setCoupling(rsTmp);
      stls.push_back(StlsCSR(inTmp));
    }
  }
  assert(stls.size() == NPOINTS);
  // Setup derivative dependency in the StlsCSR objects  
  for (size_t i = 0; i < stls.size(); ++i) {
    switch (i) {
    case RS_DOWN_THETA_DOWN: case RS_DOWN_THETA: case RS_DOWN_THETA_UP:
      stls[i].setDrsData(stls[i + 1], stls[i + 2],
			 StlsCSR::Derivative::FORWARD); break;
    case RS_THETA_DOWN: case RS_THETA: case RS_THETA_UP:
      stls[i].setDrsData(stls[i + 1], stls[i - 1],
			 StlsCSR::Derivative::CENTERED); break;
    case RS_UP_THETA_DOWN: case RS_UP_THETA: case RS_UP_THETA_UP:
      stls[i].setDrsData(stls[i - 1], stls[i - 2],
			 StlsCSR::Derivative::BACKWARD); break;
    }
  }
  for (size_t i = 0; i < stls.size(); ++i) {
    switch (i) {
    case RS_DOWN_THETA_DOWN: case RS_THETA_DOWN: case RS_UP_THETA_DOWN:
      stls[i].setDThetaData(stls[i + NRS], stls[i + 2 * NRS],
			    StlsCSR::Derivative::FORWARD); break;
    case RS_DOWN_THETA: case RS_THETA: case RS_UP_THETA:
      stls[i].setDThetaData(stls[i + NRS], stls[i - NRS],
			    StlsCSR::Derivative::CENTERED); break;
    case RS_DOWN_THETA_UP: case RS_THETA_UP: case RS_UP_THETA_UP:
      stls[i].setDThetaData(stls[i - NRS], stls[i - 2 * NRS],
			 StlsCSR::Derivative::BACKWARD); break;
    }
  }
}

// Add a public call to compute where we do the initializations
int StructProp::compute() {
  try {
    doIterations();
    computed = true;
    return 0;
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    return 1;
  }
}

void StructProp::doIterations() {
  const int maxIter = stls[0].in.getNIter();
  const double minErr = stls[0].in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  for (auto& s : stls) { s.initialGuess(); }
  // Iteration to solve for the structural properties
  while (counter < maxIter+1 && err > minErr ) {
    // Compute new solution and error
#pragma omp parallel
    {
      #pragma omp for
      for (auto& s : stls) {
	s.computeSsf();
	s.computeSlfcStls();
      }
      #pragma omp for
      for (auto& s : stls) {
	s.computeSlfc();
	s.updateSolution();
      }
    }
    counter++;
    // Compute the error only for the central state point (rs, theta)
    err = stls[RS_THETA].computeError();
  }
  printf("Alpha = %.5e, Residual error "
	 "(structural properties) = %.5e\n", stls[RS_THETA].alpha, err);
}

void StructProp::setAlpha(const double& alpha) {
  for (auto& s : stls) { s.setAlpha(alpha); }
}

const vector<double>& StructProp::getBase(function<double(const StlsCSR&)> f) const {
  for (size_t i = 0; i < NPOINTS; ++i) {
    outVector[i] = f(stls[i]);
  }
  return outVector; 
}

vector<double> StructProp::getCouplingParameters() const {
  return getBase([&](const StlsCSR& s){ return s.in.getCoupling();});
}

vector<double> StructProp::getDegeneracyParameters() const {
  return getBase([&](const StlsCSR& s){return s.in.getDegeneracy();});
}

vector<double> StructProp::getInternalEnergy() const  {
  return getBase([&](const StlsCSR& s){return s.getInternalEnergy();});
}

vector<double> StructProp::getFreeEnergyIntegrand() const  {
  return getBase([&](const StlsCSR& s){return s.getFreeEnergyIntegrand();});
}

// -----------------------------------------------------------------
// StlsCSR class
// -----------------------------------------------------------------

void StlsCSR::setDrsData(StlsCSR &stlsRsUp,
			 StlsCSR &stlsRsDown,
			 const Derivative &dTypeRs) {
  this->slfcStlsRsUp = &stlsRsUp.slfcStls;
  this->slfcStlsRsDown = &stlsRsDown.slfcStls;
  this->dTypeRs = dTypeRs;
}

void StlsCSR::setDThetaData(StlsCSR &stlsThetaUp,
			    StlsCSR &stlsThetaDown,
			    const Derivative &dTypeTheta) {
  this->slfcStlsThetaUp = &stlsThetaUp.slfcStls;
  this->slfcStlsThetaDown = &stlsThetaDown.slfcStls;
  this->dTypeTheta = dTypeTheta;
}

double StlsCSR::getDerivative(const vector<double>& f,
			      const size_t& idx,
			      const Derivative& type) {
  switch(type) {
  case BACKWARD:
    assert(idx >= 2);
    return getDerivative(f[idx], f[idx - 1], f[idx - 2], type);
    break;
  case CENTERED:
    assert(idx >= 1 && idx < f.size() - 1);
    return getDerivative(f[idx], f[idx + 1], f[idx - 1], type);
    break;
  case FORWARD:
    assert(idx < f.size() - 2);
    return getDerivative(f[idx], f[idx + 1], f[idx + 2], type);
    break;
  default:
    assert(false);
    return -1;
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
  default:
    assert(false);
    return -1;
    break;
  }
}

void StlsCSR::computeSlfcStls() {
  Stls::computeSlfc();
  slfcStls = slfcNew;
}

void StlsCSR::computeSlfc() {
  // Check that alpha has been set to a value that is not the default
  assert(alpha != DEFAULT_ALPHA);
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
  const double a_dt = alpha * theta / (3.0 * dTheta);
  const size_t nx = wvg.size();
  // Wave-vector derivative
  slfcNew[0] -= a_dx * wvg[0] * getDerivative(slfcStls, 0, FORWARD);
  for (size_t i = 1; i < nx - 1; ++i) {
    slfcNew[i] -= a_dx * wvg[i] * getDerivative(slfcStls, i, CENTERED);
  }
  slfcNew[nx - 1] -= a_dx * wvg[nx - 1] * getDerivative(slfcStls, nx - 1, BACKWARD);
  // Coupling parameter contribution
  if (rs > 0.0) {
    for (size_t i = 0; i < nx; ++i) {
      slfcNew[i] -= a_drs * getDerivative(slfcStls[i], rsUp[i], rsDown[i], dTypeRs);
    }
  }
  // Degeneracy parameter contribution
  if (theta > 0.0) {
    for (size_t i = 0; i < nx; ++i) {
      slfcNew[i] -= a_dt * getDerivative(slfcStls[i], thetaUp[i], thetaDown[i], dTypeTheta);
    }
  }
}

double StlsCSR::getInternalEnergy() const {
  return computeInternalEnergy(wvg, ssf, in.getCoupling());
}

double StlsCSR::getFreeEnergyIntegrand() const {
  return computeInternalEnergy(wvg, ssf, 1.0);
}
