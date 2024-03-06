#include "util.hpp"
#include "numerics.hpp"
#include "input.hpp"
#include "qvs.hpp"

using namespace std;
using namespace vecUtil;

// -----------------------------------------------------------------
// qVSStls class
// -----------------------------------------------------------------

double qVSStls::computeAlpha() {
  //  Compute the free energy integrand
  thermoProp.compute<qVSStls>(in);
  // Free energy
  const vector<double> freeEnergyData = thermoProp.getFreeEnergyData();
  const double& fxcr = freeEnergyData[1];
  const double& fxcrr = freeEnergyData[2];
  const double& fxct = freeEnergyData[3];
  const double& fxctt = freeEnergyData[4];
  const double& fxcrt = freeEnergyData[5];
  // Q
  const vector<double> QData = thermoProp.getQData();
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
  const auto& qstls = thermoProp.getStructProp<qStlsCSR>();
  adr = qstls.getAdr();
  ssf = qstls.getSsf();
}

// -----------------------------------------------------------------
// qThermoProp class
// -----------------------------------------------------------------

vector<double> qThermoProp::getQData() const {
  // COMMENTS:
  // 1) We don't need to define a qStructProp object here, structProp already
  // has a method called getQ because it is of type qStructProp. This might
  // sound a bit confusing, but is an effect of using the templates: qThermoProp
  // is equivalent to ThermoPropBase with StructProp = qStructProp,
  // just look at line 124 in qvs.hpp. If you feel that you are getting lost
  // try to attach a debugger and check the types of the variables inside
  // qThermoProp.
  // 2) QAdder diverges for rs = 0, this can cause you problems when computing
  // the derivative with respect to the coupling parameter. To circumvent this
  // problem (which appears also with the internal energy), check what I have
  // done in ThermoPropBase::getInternalEnergyData(). Essentially, instead of
  // using the derivatives of QAdder, work with the derivative of rs * QAdder,
  // which is finite. My suggestion is to remove the 1/rs term from the definition
  // of QAdder and then re-introduce it only here when you need q and qt, e.g.
  // replace q = qVec[SIdx::RS_THETA] with q = qVec[SIdx::RS_THETA] / rs[SIdx::RS_THETA]
  // Note that q and qt are not affected by the problem of the divergence because
  // we never compute alpha directly for rs = 0, we always start at rs = drs.
  // QAdder
  const std::vector<double> qVec = structProp.getQ();
  const double q = qVec[SIdx::RS_THETA];
  // QAdder derivative with respect to the coupling parameter
  double qr;
  {
    const std::vector<double> rs = structProp.getCouplingParameters(); 
    const double drs = rs[SIdx::RS_UP_THETA] - rs[SIdx::RS_THETA];
    const double& q0 = qVec[SIdx::RS_UP_THETA];
    const double& q1 = qVec[SIdx::RS_DOWN_THETA];
    qr = rs[SIdx::RS_THETA] * (q0 - q1) / (2.0 * drs);
  }
  // QAdder derivative with respect to the degeneracy parameter
  double qt;
  {
    const std::vector<double> theta = structProp.getDegeneracyParameters();
    const double dt = theta[SIdx::RS_THETA_UP] - theta[SIdx::RS_THETA];
    const double q0 = qVec[SIdx::RS_THETA_UP];
    const double q1 = qVec[SIdx::RS_THETA_DOWN];
    qt = theta[SIdx::RS_THETA] * (q0 - q1) / (2.0 * dt);
  }
  return vector<double>({q, qr, qt});
}

// -----------------------------------------------------------------
// qStructProp class
// -----------------------------------------------------------------

void qStructProp::doIterations() {
  const auto& in = csr[0].in;
  const int maxIter = in.getNIter();
  const double minErr = in.getErrMin();
  double err = 1.0;
  int counter = 0;
  // Define initial guess
  for (auto& c : csr) { c.initialGuess(); }
  // Iteration to solve for the structural properties
  while (counter < maxIter+1 && err > minErr ) {
    for (auto& c : csr) {
      c.computeAdrStls();
      c.computeSsf();
      c.computeAdr();
      c.updateSolution();
    }
    counter++;
    // Compute the error only for the central state point (rs, theta)
    err = csr[RS_THETA].computeError();
  }
  printf("Alpha = %.5e, Residual error "
	 "(structural properties) = %.5e\n", csr[RS_THETA].alpha, err);
}

vector<double> qStructProp::getQ() const  {
  return getBase([&](const qStlsCSR& q) -> double {
    return q.getQAdder();
  });
}

// -----------------------------------------------------------------
// qStlsCSR class
// -----------------------------------------------------------------

void qStlsCSR::computeAdrStls() {
  Qstls::computeAdr();
  lfc = adr;
}

double qStlsCSR::getDerivative(const Vector2D& f,
			       const int &l,
			       const size_t& idx,
			       const Derivative& type) {
  switch(type) {
  case BACKWARD:
    assert(idx >= 2);
    return CSR::getDerivative(f(idx,l), f(idx - 1,l), f(idx - 2,l), type);
    break;
  case CENTERED:
    assert(idx >= 1 && idx < f.size() - 1);
    return CSR::getDerivative(f(idx,l), f(idx + 1,l), f(idx - 1,l), type);
    break;
  case FORWARD:
    assert(idx < f.size() - 2);
    return CSR::getDerivative(f(idx,l), f(idx + 1,l), f(idx + 2,l), type);
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
  for (int l = 0; l < nl; ++l) {
    // Wave-vector derivative contribution
    adr(0,l) -= a_dx * wvg[0] * getDerivative(lfc, l, 0, FORWARD);
    for (size_t i = 1; i < nx - 1; ++i) {
      adr(i,l) -= a_dx * wvg[i] * getDerivative(lfc, l, i, CENTERED);
    }
    adr(nx - 1,l) -= a_dx * wvg[nx - 1] * getDerivative(lfc, l, nx - 1, BACKWARD);
    // Coupling parameter contribution
    if (rs > 0.0) {
      for (size_t i = 0; i < nx; ++i) {
        adr(i,l) -= a_drs * CSR::getDerivative(lfc(i,l), rsUp(i,l), rsDown(i,l), dTypeRs);
      }
    }
    // Degeneracy parameter contribution
    if (theta > 0.0) {
      for (size_t i = 0; i < nx; ++i) {
        adr(i,l) -= a_dt * CSR::getDerivative(lfc(i,l), thetaUp(i,l), thetaDown(i,l), dTypeTheta);
      }
    }
    // Extra 1/3 term present in the new adr for qVS
    for (size_t i = 0; i < nx; ++i) {
      adr(i,l) += alpha/3.0 * lfc(i,l);
    }
  }
}

double qStlsCSR::getQAdder() const {
  Integrator1D itg1(in.getIntError());
  Integrator2D itg2(in.getIntError());
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  const Interpolator1D ssfItp(wvg, ssf);
  QAdder QTmp(in.getCoupling(), in.getDegeneracy(), mu, wvg, 
              itgGrid, itg1, itg2, ssfItp); 
  return QTmp.get();
}

// // -----------------------------------------------------------------
// // QAdder class
// // -----------------------------------------------------------------

// SSF interpolation
double QAdder::ssf(const double& y) const {
  return interp.eval(y);
}

// Denominator integrand
double QAdder::integrandDenominator(const double y) const {
  const double y2 = y*y;
  return 1.0/(exp(y2/Theta - mu) + 1.0);
}

// Numerator integrand1
double QAdder::integrandNumerator1(const double q) const {
  const double w = itg2.getX();
  if (w == 0 || w == 2*q) { return 0; };
  const double w2 = w*w;
  const double w3 = w2*w;
  const double logarg = (w + 2*q)/(w - 2*q);
  return q/(exp(q*q/Theta - mu) + 1.0) * q/w3 * (q/w*log(logarg) - 1.0);
}

// Numerator integrand2
double QAdder::integrandNumerator2(const double w) const {
  if (w == 0.0) return 0.0;
  return w * (ssf(w) - 1.0);
}

// Denominator integral
void QAdder::getIntDenominator(double &res) const {
  // COMMENT: back() refers to the last element in the vector, front() to
  // the first one. You seem to use them backwards when defining yMin
  // and yMax
  auto yMin = wvg.back();
  auto yMax = wvg.front();
  auto func = [&](double y)->double{return integrandDenominator(y);};
  itg1.compute(func, yMin, yMax);
  res = itg1.getSolution();
}

// Get at finite temperature
double QAdder::get() const {
  double Denominator;
  getIntDenominator(Denominator);
  // COMMENT: back() refers to the last element in the vector, front() to
  // the first one. You seem to use them backwards when defining yMin
  // and yMax
  auto wMin = wvg.back();
  auto wMax = wvg.front();
  auto yMin = wMin;
  auto yMax = wMax;
  auto func1 = [&](const double& q)->double{return integrandNumerator1(q);};
  auto func2 = [&](const double& w)->double{return integrandNumerator2(w);};
  itg2.compute(func1, func2, yMin, yMax, yMax, wMax, itgGrid);
  return 12.0 / (M_PI * lambda * rs) * itg2.getSolution()/Denominator;
}
