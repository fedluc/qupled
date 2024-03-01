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

double qVSStls::alphaDifference(const double& alphaTmp) {
  alpha = alphaTmp;
  qthermoProp.setAlpha(alpha);
  const double alphaTheoretical = computeAlpha();
  return alpha - alphaTheoretical;
}

double qVSStls::computeAlpha() {
  // Compute the free energy integrand
  qthermoProp.compute(in);
  // Free energy
  const vector<double> freeEnergyData = qthermoProp.getFreeEnergyData();
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
  const auto& Adr = qthermoProp.getStructProp();
  adr = Adr.getAdr();
  ssf = Adr.getSsf();
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
      fxciTmp = qstructProp.getFreeEnergyIntegrand();
    }
    else if (rs < in.getCoupling()) {
      if (rs == 0.0 || fxcIntegrand[THETA][i] != Inf) { continue; }
      printf("Free energy integrand calculation, solving VS-STLS scheme for rs = %.5f:\n", rs);
      inTmp.setCoupling(rs);
      qVSStls vsstlsTmp(inTmp, *this);
      vsstlsTmp.compute();
      fxciTmp = vsstlsTmp.getThermoProp().qstructProp.getFreeEnergyIntegrand();
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

const qStlsCSR& qThermoProp::getStructProp() {
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
  // // Q
  // const double q = qstructProp.getQ()[SIdx::RS_THETA];
  // // Q derivative with respect to the coupling parameter
  // double qr;
  // {
  //   const vector<double> rs = qstructProp.getCouplingParameters(); 
  //   const double drs = rs[SIdx::RS_UP_THETA] - rs[SIdx::RS_THETA];
  //   const double& q0 = qstructProp.getQ()[SIdx::RS_UP_THETA];
  //   const double& q1 = qstructProp.getQ()[SIdx::RS_DOWN_THETA];
  //   qr = (q0 - q1) / (2.0 * drs);
  // }
  // // Q derivative with respect to the degeneracy parameter
  // double qt;
  // {
  //   const vector<double> theta = qstructProp.getDegeneracyParameters();
  //   const double dt = theta[SIdx::RS_THETA_UP] - theta[SIdx::RS_THETA];
  //   const double q0 = qstructProp.getQ()[SIdx::RS_THETA_UP];
  //   const double q1 = qstructProp.getQ()[SIdx::RS_THETA_DOWN];
  //   qt = theta[SIdx::RS_THETA] * (q0 - q1) / (2.0 * dt);
  // }
  // return vector<double>({q, qr, qt});
  // PLACEHOLDER VALUE, REMOVE WHEN YOU HAVE A WORKING IMPLEMENTATION FOR THE Q TERM
  return vector<double>();
}



// -----------------------------------------------------------------
// qStructProp class
// -----------------------------------------------------------------

qStructProp::qStructProp(const VSStlsInput &in) : StructProp(in),
						  adrIsInitialized(false) {
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
      Adr.push_back(qStlsCSR(inTmp)); 
    }
  }
  assert(Adr.size() == NPOINTS);
  // Setup derivative dependency in the StlsCSR objects  
  for (size_t i = 0; i < Adr.size(); ++i) {
    switch (i) {
    case RS_DOWN_THETA_DOWN: case RS_DOWN_THETA: case RS_DOWN_THETA_UP:
      Adr[i].setDrsData(Adr[i + 1], Adr[i + 2],
			qStlsCSR::Derivative::FORWARD); break;
    case RS_THETA_DOWN: case RS_THETA: case RS_THETA_UP:
      Adr[i].setDrsData(Adr[i + 1], Adr[i - 1],
			qStlsCSR::Derivative::CENTERED); break;
    case RS_UP_THETA_DOWN: case RS_UP_THETA: case RS_UP_THETA_UP:
      Adr[i].setDrsData(Adr[i - 1], Adr[i - 2],
			 qStlsCSR::Derivative::BACKWARD); break;
    }
  }
  for (size_t i = 0; i < Adr.size(); ++i) {
    switch (i) {
    case RS_DOWN_THETA_DOWN: case RS_THETA_DOWN: case RS_UP_THETA_DOWN:
      Adr[i].setDThetaData(Adr[i + NRS], Adr[i + 2 * NRS],
			   qStlsCSR::Derivative::FORWARD); break;
    case RS_DOWN_THETA: case RS_THETA: case RS_UP_THETA:
      Adr[i].setDThetaData(Adr[i + NRS], Adr[i - NRS],
			   qStlsCSR::Derivative::CENTERED); break;
    case RS_DOWN_THETA_UP: case RS_THETA_UP: case RS_UP_THETA_UP:
      Adr[i].setDThetaData(Adr[i - NRS], Adr[i - 2 * NRS],
			   qStlsCSR::Derivative::BACKWARD); break;
    }
  }
}

int qStructProp::compute() {
  try {
    if (!adrIsInitialized) {
      for (auto& q : Adr) { q.init(); }
      adrIsInitialized = true;
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
  const auto& in = Adr[0].in;
  const int maxIter = in.getNIter();
  const int ompThreads = in.getNThreads();
  const double minErr = in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  for (auto& q : Adr) { q.initialGuess(); }
  // Iteration to solve for the structural properties
  const bool useOMP = ompThreads > 1;
  while (counter < maxIter+1 && err > minErr ) {
    // Compute new solution and error
    #pragma omp parallel num_threads(ompThreads) if (useOMP)
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

const vector<double>& qStructProp::getBase(std::function<double(const qStlsCSR&)> f) const {
  for (size_t i = 0; i < NPOINTS; ++i) {
    outVector[i] = f(Adr[i]);
  }
  return outVector; 
}

vector<double> qStructProp::getCouplingParameters() const {
  return getBase([&](const qStlsCSR& q) -> double {
    return q.in.getCoupling();
  });
}

vector<double> qStructProp::getDegeneracyParameters() const {
  return getBase([&](const qStlsCSR& q) -> double {
    return q.in.getDegeneracy();
  });
}

vector<double> qStructProp::getInternalEnergy() const  {
  return getBase([&](const qStlsCSR& q) -> double {
    return q.getInternalEnergy();
  });
}

vector<double> qStructProp::getFreeEnergyIntegrand() const  {
  return getBase([&](const qStlsCSR& q) -> double {
    return q.getFreeEnergyIntegrand();
  });
}

vector<double> qStructProp::getQ() const  {
  return getBase([&](const qStlsCSR& q) -> double {
    return q.getQ();
  });
}

// -----------------------------------------------------------------
// qStlsCSR class
// -----------------------------------------------------------------

QstlsInput qStlsCSR::VStoQStlsInput(const VSStlsInput& in) const {
  StlsInput tmp = in;
  return static_cast<QstlsInput>(tmp);
}


void qStlsCSR::computeAdrStls() {
  Qstls::computeAdr();
  lfc = adr;
}

double qStlsCSR::getqDerivative(const Vector2D& f,
           const int &l,
			     const size_t& idx,
			     const Derivative& type) {
  // NOTE: If T does not have an operator[] this method would not compile
  switch(type) {
  case BACKWARD:
    assert(idx >= 2);
    return getDerivative(f(idx,l), f(idx - 1,l), f(idx - 2,l), type);
    break;
  case CENTERED:
    assert(idx >= 1 && idx < f.size() - 1);
    return getDerivative(f(idx,l), f(idx + 1,l), f(idx - 1,l), type);
    break;
  case FORWARD:
    assert(idx < f.size() - 2);
    return getDerivative(f(idx,l), f(idx + 1,l), f(idx + 2,l), type);
    break;
  default:
    assert(false);
    return -1;
    break;
  }
}

void qStlsCSR::computeAdr() {
  // Check that alpha has been set to a value that is not the default
  assert(alpha != DEFAULT_ALPHA);
  // Derivative contributions
  const double& rs = in.getCoupling();
  //const double& theta = in.getDegeneracy();
  const double& theta = 0.0;
  const double& dx = in.getWaveVectorGridRes();
  const double& drs = in.getCouplingResolution();
  const double& dTheta = in.getDegeneracyResolution();
  const Vector2D& rsUp = *lfcRsUp;
  const Vector2D& rsDown = *lfcRsDown;
  const Vector2D& thetaUp = *lfcThetaUp;
  const Vector2D& thetaDown = *lfcThetaDown;
  const double a_drs = alpha * rs / (6.0 * drs);
  const double a_dx = alpha/(6.0 * dx);
  const double a_dt = alpha * theta / (3.0 * dTheta);
  const size_t nx = wvg.size();
  const int nl = adrFixed.size(1);
  // Wave-vector derivative contribution
  for (int l = 0; l < nl; ++l) {
    adr(0,l) -= a_dx * wvg[0] * getqDerivative(lfc, l, 0, FORWARD);
    for (size_t i = 1; i < nx - 1; ++i) {
      adr(i,l) -= a_dx * wvg[i] * getqDerivative(lfc, l, i, CENTERED);
    }
    adr(nx - 1,l) -= a_dx * wvg[nx - 1] * getqDerivative(lfc, l, nx - 1, BACKWARD);
    // Coupling parameter contribution
    if (rs > 0.0) {
      for (size_t i = 0; i < nx; ++i) {
        adr(i,l) -= a_drs * getDerivative(lfc(i,l), rsUp(i,l), rsDown(i,l), dTypeRs);
      }
    }
    // Degeneracy parameter contribution
    if (theta > 0.0) {
      for (size_t i = 0; i < nx; ++i) {
        adr(i,l) -= a_dt * getDerivative(lfc(i,l), thetaUp(i,l), thetaDown(i,l), dTypeTheta);
      }
    }
    // Extra 1/3 term present in the new adr for qVS
    for (size_t i = 0; i < nx; ++i) {
      adr(i,l) += alpha/3.0 * lfc(i,l);
    }
  }
}

// THIS PART STILL NEEDS WORK
double qStlsCSR::getQ() const {
  return QInstance.computeQ(wvg, ssf, in.getCoupling(), in.getDegeneracy());
}

// Integrands for the fixed component
double Q::integrandDenominator(const double y) const {
  const double y2 = y*y;
  return 1.0/(exp(y2/Theta - mu) + 1.0);
}

double Q::integrandNumerator(const double q,
			    const double w) const {
  const double Theta = in.getDegeneracy();
  //const double q = itg.getX();
  if (q == 0 || w == 0) { return 0; };
  const double w2 = w*w;
  const double w3 = w2*w;
  const double logarg = (w + 2*q)/(w - 2*q);
  return q/(exp(q*q/Theta - mu) + 1.0) * q/w3 * (q/w*log(logarg) - 1.0);
}

// Denominator integral
void Q::getIntDenominator(const vector<double> wvg,
        double &res) const {
  auto yMin = wvg.front();
  auto yMax = wvg.back();
  auto func = [&](double y)->double{return integrandDenominator(y);};
  itg.compute(func, yMin, yMax);
  res = itg.getSolution();
}
// Numerator integral without ssf
void Q::getIntNumerator(const vector<double> wvg,
        std::vector<double> &res) const{
  auto qMin = wvg.front();
  auto qMax = wvg.back();
  res.resize(wvg.size());
  for (int i = 0; i < wvg.size(); ++i) {
    double w = wvg[i];
    auto func = [&](double q)->double{return integrandNumerator(q,w);};
    itg.compute(func, qMin, qMax);
    res[i] = itg.getSolution();
  }
}

void Q::getTotalNumerator(const vector<double> wvg,
        std::vector<double> &NumPart,
        std::vector<double> ssf,
        double &Numerator) const{
  getIntNumerator(wvg, NumPart);
  double wMin = wvg.front();
  double wMax = wvg.back();
  assert(wvg.size() == NumPart.size() && wvg.size() == ssf.size());
  interp.reset(wvg[0], NumPart[0], wvg.size());

  auto totalIntegrand = [this, &wvg, &ssf](double w) -> double {
    double interpolatedValue = this->interp.eval(w);
    auto it = std::find(wvg.begin(), wvg.end(), w);
    if (it == wvg.end()) {
      return 0.0; 
    }
    size_t index = std::distance(wvg.begin(), it);
    double WssfMinusOne = w * (ssf[index] - 1);
    return WssfMinusOne * interpolatedValue;
    };

  itg.compute(totalIntegrand, wMin, wMax);
  Numerator = itg.getSolution();
}

double Q::computeQ(const vector<double> &wvg, 
                  const std::vector<double> ssf, 
                  const double rs, 
                  const double Theta) {
  getIntDenominator(wvg, Denominator);
return 12.0 / (M_PI * lambda * rs) * (Numerator / Denominator);
}

