#include <omp.h>
#include "util.hpp"
#include "numerics.hpp"
#include "input.hpp"
#include "qvs.hpp"


using namespace std;
using namespace numUtil;
using namespace vecUtil;
using namespace thermoUtil;
using namespace parallelUtil;

// -----------------------------------------------------------------
// qVSStls class
// -----------------------------------------------------------------

int qVSStls::compute() {
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

// qStls iterations
void qVSStls::doIterations() {
  auto func = [this](double alphaTmp)->double{return alphaDifference(alphaTmp);};
  SecantSolver rsol(in.getErrMinAlpha(), in.getNIterAlpha());
  rsol.solve(func, in.getAlphaGuess());
  if (!rsol.success()) {
    MPI::throwError("VSStls: the root solver did not converge to the desired accuracy.");
  }
  alpha = rsol.getSolution();
  if (verbose) { cout << "Free parameter = " << alpha << endl; }
  updateSolution();
}

double qVSStls::alphaDifference(const double& alphaTmp) {
  alpha = alphaTmp;
  thermoProp.setAlpha(alpha);
  const double alphaTheoretical = computeAlpha();
  return alpha - alphaTheoretical;
}

double qVSStls::computeAlpha() {
  // Compute the free energy integrand
  qthermoProp.compute(in);
  // Free energy
  const vector<double> freeEnergyData = thermoProp.getFreeEnergyData();
  const double& fxc = freeEnergyData[0];
  const double& fxcr = freeEnergyData[1];
  const double& fxcrr = freeEnergyData[2];
  const double& fxct = freeEnergyData[3];
  const double& fxctt = freeEnergyData[4];
  const double& fxcrt = freeEnergyData[5];
  // Q
  const vector<double> QData = qthermoProp.getQData();
  const double& Q = QData[0];
  const double& Qr = QData[1];
  const double& Qt = QData[2];
  // Alpha
  double numer = Q - (1.0/6.0) * fxcrr + (1.0/3.0) * fxcr;
  double denom =  Q + (1.0/3.0) * Qr;
  if (in.getDegeneracy() > 0.0) {
    numer += - (2.0/3.0) * fxctt
             - (2.0/3.0) * fxcrt
             + (1.0/3.0) * fxct;
    denom += (2.0/3.0) * Qt;
  }
  return numer/denom;
}

void qVSStls::updateSolution() {
  // Update the structural properties used for output
  const auto& Adr = qthermoProp.getStructProp<std::vector<double>>();
  adr = Adr.getAdr();
  ssf = Adr.get();
}

// -----------------------------------------------------------------
// qThermoProp class
// -----------------------------------------------------------------

void qThermoProp::setAlpha(const double& alpha) {
  qstructProp.setAlpha(alpha);
}

void qThermoProp::compute(const VSStlsInput& in) {
  // Recursive calls to solve the VS-STLS scheme for all state points
  // with coupling parameter smaller than rs
  const double nrs = rsGrid.size();
  VSStlsInput inTmp = in;
  vector<double> fxciTmp(StructProp::NPOINTS);
  for (size_t i = 0; i < nrs; ++i) {
    const double& rs = rsGrid[i];
    if (equalTol(rs, in.getCoupling())) {
      qstructProp.compute();
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

template<typename T = std::vector<double>>
const qStlsCSR<T>& qThermoProp::getStructProp() {
  if (!qstructProp.isComputed()) { qstructProp.compute(); }
  if (isZeroCoupling && isZeroDegeneracy) {
     return qstructProp.getAdrStls(SIdx::RS_DOWN_THETA_DOWN); 
  }
  if (!isZeroCoupling && isZeroDegeneracy) {
     return qstructProp.getAdrStls(SIdx::RS_THETA_DOWN); 
  }
  if (isZeroCoupling && !isZeroDegeneracy) {
     return qstructProp.getAdrStls(SIdx::RS_DOWN_THETA); 
  }
  return qstructProp.getAdrStls(SIdx::RS_THETA); 
}


vector<double> qThermoProp::getQData() const {
  // Q
  const double q = qstructProp.getQ()[SIdx::RS_THETA];
  // Q derivative with respect to the coupling parameter
  double qr;
  {
    const vector<double> rs = structProp.getCouplingParameters(); 
    const double drs = rs[SIdx::RS_UP_THETA] - rs[SIdx::RS_THETA];
    const double& q0 = qstructProp.getQ()[SIdx::RS_UP_THETA];
    const double& q1 = qstructProp.getQ()[SIdx::RS_DOWN_THETA];
    qr = (q0 - q1) / (2.0 * drs);
  }
  // Q derivative with respect to the degeneracy parameter
  double qt;
  {
    const vector<double> theta = structProp.getDegeneracyParameters();
    const double dt = theta[SIdx::RS_THETA_UP] - theta[SIdx::RS_THETA];
    const double q0 = qstructProp.getQ()[SIdx::RS_THETA_UP];
    const double q1 = qstructProp.getQ()[SIdx::RS_THETA_DOWN];
    qt = theta[SIdx::RS_THETA] * (q0 - q1) / (2.0 * dt);
  }
  return vector<double>({q, qr, qt});
}



// -----------------------------------------------------------------
// qStructProp class
// -----------------------------------------------------------------

StructProp::StructProp(const VSStlsInput &in) : stlsIsInitialized(false),
						computed(false),
						outVector(NPOINTS) {
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

int qStructProp::compute() {
  try {
    if (!stlsAdrIsInitialized) {
      for (auto& q : Adr) { q.init(); }
      stlsAdrIsInitialized = true;
    }
    doIterations();
    computed = true;
    return 0;
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    return 1;
  }
}

void qStructProp::doIterations() {
  const int maxIter = Adr[0].in.getNIter();
  const double minErr = Adr[0].in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  for (auto& q : Adr) { q.initialGuess(); }
  // Iteration to solve for the structural properties
  while (counter < maxIter+1 && err > minErr ) {
    // Compute new solution and error
#pragma omp parallel
    {
      #pragma omp for
      for (auto& q : Adr) {
  q.computeAdrStls();     
	q.computeSsf();
      }
      #pragma omp for
      for (auto& q : Adr) {
	q.computeAdr();
	q.updateSolution();
      }
    }
    counter++;
    // Compute the error only for the central state point (rs, theta)
    err = Adr[RS_THETA].computeError();
  }
  printf("Alpha = %.5e, Residual error "
	 "(structural properties) = %.5e\n", Adr[RS_THETA].alpha, err);
}

void qStructProp::setAlpha(const double& alpha) {
  for (auto& q : Adr) { q.setAlpha(alpha); }
}

const vector<double>& qStructProp::getBase(std::function<double(const qStlsCSR<std::vector<double>>&)> f) const {
  for (size_t i = 0; i < NPOINTS; ++i) {
    outVector[i] = f(Adr[i]);
  }
  return outVector; 
}

vector<double> qStructProp::getCouplingParameters() const {
  return getBase([&](const qStlsCSR<std::vector<double>>& q) -> double {
    return q.in.getCoupling();
  });
}

vector<double> qStructProp::getDegeneracyParameters() const {
  return getBase([&](const qStlsCSR<std::vector<double>>& q) -> double {
    return q.in.getDegeneracy();
  });
}

vector<double> qStructProp::getInternalEnergy() const  {
  return getBase([&](const qStlsCSR<std::vector<double>>& q) -> double {
    return q.getInternalEnergy();
  });
}

vector<double> qStructProp::getFreeEnergyIntegrand() const  {
  return getBase([&](const qStlsCSR<std::vector<double>>& q) -> double {
    return q.getFreeEnergyIntegrand();
  });
}

vector<double> qStructProp::getQ() const  {
  return getBase([&](const qStlsCSR<std::vector<double>>& q) -> double {
    return q.getQ();
  });
}

// -----------------------------------------------------------------
// qCSR class
// -----------------------------------------------------------------

template<typename T>
void CSR<T>::setDrsData(CSR<T> &csrRsUp,
			CSR<T> &csrRsDown,
			const Derivative &dTypeRs) {
  this->lfcRsUp = &csrRsUp.lfc;
  this->lfcRsDown = &csrRsDown.lfc;
  this->dTypeRs = dTypeRs;
}

template<typename T>
void CSR<T>::setDThetaData(CSR<T> &csrThetaUp,
			   CSR<T> &csrThetaDown,
			   const Derivative &dTypeTheta) {
  this->lfcThetaUp = &csrThetaUp.lfc;
  this->lfcThetaDown = &csrThetaDown.lfc;
  this->dTypeTheta = dTypeTheta;
}

template<typename T>
double CSR<T>::getDerivative(const T& f,
			     const size_t& idx,
			     const Derivative& type) {
  // NOTE: If T does not have an operator[] this method would not compile
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

template<typename T>
double CSR<T>::getDerivative(const double& f0,
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

// -----------------------------------------------------------------
// qStlsCSR class
// -----------------------------------------------------------------

template<typename T>
void qStlsCSR<T>::computeAdrStls() {
  Qstls::computeAdr();
  lfc = adr;
}

template<typename T>
void qStlsCSR<T>::computeAdr() {
  // Check that alpha has been set to a value that is not the default
  assert(alpha != DEFAULT_ALPHA);
  // Derivative contributions
  const double& rs = in.getCoupling();
  //const double& theta = in.getDegeneracy();
  const double& theta = 0.0;
  const double& dx = in.getWaveVectorGridRes();
  const double& drs = in.getCouplingResolution();
  const double& dTheta = in.getDegeneracyResolution();
  const vector<double>& rsUp = *lfcRsUp;
  const vector<double>& rsDown = *lfcRsDown;
  const vector<double>& thetaUp = *lfcThetaUp;
  const vector<double>& thetaDown = *lfcThetaDown;
  const double a_drs = alpha * rs / (6.0 * drs);
  const double a_dx = alpha/(6.0 * dx);
  const double a_dt = alpha * theta / (3.0 * dTheta);
  const size_t nx = wvg.size();
  // Wave-vector derivative
  adr[0] -= a_dx * wvg[0] * getDerivative(lfc, 0, FORWARD);
  for (size_t i = 1; i < nx - 1; ++i) {
    adr[i] -= a_dx * wvg[i] * getDerivative(lfc, i, CENTERED);
  }
  adr[nx - 1] -= a_dx * wvg[nx - 1] * getDerivative(lfc, nx - 1, BACKWARD);
  // Coupling parameter contribution
  if (rs > 0.0) {
    for (size_t i = 0; i < nx; ++i) {
      adr[i] -= a_drs * getDerivative(lfc[i], rsUp[i], rsDown[i], dTypeRs);
    }
  }
  // Degeneracy parameter contribution
  if (theta > 0.0) {
    for (size_t i = 0; i < nx; ++i) {
      adr[i] -= a_dt * getDerivative(lfc[i], thetaUp[i], thetaDown[i], dTypeTheta);
    }
  }
}

// THIS PART STILL NEEDS WORK
template<typename T>
double qStlsCSR<T>::getQ() const {

  return computeQ(wvg, ssf, in.getCoupling(), in.getDegeneracy());
}

// Get fixed component
void Q::get(vector<double> &wvg,
		   Vector3D &res) const {
  const int nx = wvg.size();
  const int nl = res.size(1);
  if ( x == 0.0 ) { res.fill(0.0); };
  const double x2 = x*x;
  auto it = find(wvg.begin(), wvg.end(), x);
  assert(it != wvg.end());
  size_t ix = distance(wvg.begin(), it);
  for (int l = 0; l < nl; ++l){
    for (int i = 0; i < nx; ++i) {
      const double xq = x*wvg[i];
      auto tMin = [&]()->double{return x2 - xq;};
      auto tMax = [&]()->double{return x2 + xq;};
      auto func1 = [&](double q)->double{return integrandnum(q, l);};
      auto func2 = [&](double t)->double{return integranddenom(t, wvg[i], l);};
      itg.compute(func1, func2, qMin, qMax, tMin, tMax, itgGrid);
      res(ix, l, i) = itg.getSolution();
    }
  }
}

// Integrands for the fixed component
double Q::integranddenom(const double q,
			    const double l) const {
  if (l == 0) return q/(exp(q*q/Theta - mu) + exp(-q*q/Theta + mu) + 2.0);
  return 1.0/(exp(q*q/Theta - mu) + 1.0);
}

double Q::integrandnum(const double t,
			    const double y,
			    const double l) const {
  const double q = itg.getX();
  if (q == 0 || t == 0 || y == 0) { return 0; };
  const double x2 = x*x;
  const double y2 = y*y;
  const double y3 = y2*y;
  const double logarg = (y + 2*q)/(y - 2*q);
  return q/(exp(q*q/Theta - mu) + 1.0) * q/y3 * (q/y*log(logarg) - 1.0);
}