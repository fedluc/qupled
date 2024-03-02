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
  parallelUtil::MPI::throwError("qThermoProp::getQData() temporarily disabled"); // REMOVE WHEN YOU HAVE A WORKING IMPLEMENTATION
  return vector<double>();
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
    return q.getQ();
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

double qStlsCSR::getQ() const {
  //return QInstance.computeQ(wvg, ssf, in.getCoupling(), in.getDegeneracy());
  // COMMENT: Instead of calling Qinstance.computeQ simply create a local
  // object of type Q and call computeQ on it. For example:
  // Q q(<constructor parameters>);
  // return q.compute(...);
  parallelUtil::MPI::throwError("qStlsCSR::getQ() temporarily disabled"); // REMOVE WHEN YOU HAVE A WORKING IMPLEMENTATION
  return numUtil::Inf;
}

// // -----------------------------------------------------------------
// // Q class
// // -----------------------------------------------------------------

// // Integrands for the fixed component
// double Q::integrandDenominator(const double y) const {
//   const double y2 = y*y;
//   return 1.0/(exp(y2/Theta - mu) + 1.0);
// }

// double Q::integrandNumerator(const double q,
// 			    const double w) const {
//   const double Theta = in.getDegeneracy();
//   if (q == 0 || w == 0) { return 0; };
//   const double w2 = w*w;
//   const double w3 = w2*w;
//   const double logarg = (w + 2*q)/(w - 2*q);
//   return q/(exp(q*q/Theta - mu) + 1.0) * q/w3 * (q/w*log(logarg) - 1.0);
// }

// // Denominator integral
// void Q::getIntDenominator(const vector<double> wvg,
//         double &res) const {
//   auto yMin = wvg.front();
//   auto yMax = wvg.back();
//   auto func = [&](double y)->double{return integrandDenominator(y);};
//   itg.compute(func, yMin, yMax);
//   res = itg.getSolution();
// }
// // Numerator integral without ssf
// void Q::getIntNumerator(const vector<double> wvg,
//         std::vector<double> &res) const{
//   auto qMin = wvg.front();
//   auto qMax = wvg.back();
//   res.resize(wvg.size());
//   for (int i = 0; i < wvg.size(); ++i) {
//     double w = wvg[i];
//     auto func = [&](double q)->double{return integrandNumerator(q,w);};
//     itg.compute(func, qMin, qMax);
//     res[i] = itg.getSolution();
//   }
// }

// void Q::getTotalNumerator(const vector<double> wvg,
//         std::vector<double> &NumPart,
//         std::vector<double> ssf,
//         double &Numerator) const{
//   getIntNumerator(wvg, NumPart);
//   double wMin = wvg.front();
//   double wMax = wvg.back();
//   assert(wvg.size() == NumPart.size() && wvg.size() == ssf.size());
//   interp.reset(wvg[0], NumPart[0], wvg.size());

//   auto totalIntegrand = [this, &wvg, &ssf](double w) -> double {
//     double interpolatedValue = this->interp.eval(w);
//     auto it = std::find(wvg.begin(), wvg.end(), w);
//     if (it == wvg.end()) {
//       return 0.0; 
//     }
//     size_t index = std::distance(wvg.begin(), it);
//     double WssfMinusOne = w * (ssf[index] - 1);
//     return WssfMinusOne * interpolatedValue;
//     };

//   itg.compute(totalIntegrand, wMin, wMax);
//   Numerator = itg.getSolution();
// }

// double Q::computeQ(const vector<double> &wvg, 
//                   const std::vector<double> ssf, 
//                   const double rs, 
//                   const double Theta) {
//   getIntDenominator(wvg, Denominator);
// return 12.0 / (M_PI * lambda * rs) * (Numerator / Denominator);
// }

