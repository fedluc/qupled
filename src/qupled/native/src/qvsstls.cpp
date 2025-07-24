#include "qvsstls.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"
#include <filesystem>

using namespace std;
using namespace vecUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;
using ItgType = Integrator1D::Type;

// -----------------------------------------------------------------
// QVSStls class
// -----------------------------------------------------------------

QVSStls::QVSStls(const std::shared_ptr<const QVSStlsInput> &in_)
    : VSBase(),
      Qstls(in_, false),
      thermoProp(make_shared<QThermoProp>(in_)) {
  VSBase::thermoProp = thermoProp;
}

double QVSStls::computeAlpha() {
  //  Compute the free energy integrand
  thermoProp->compute();
  // Free energy
  const vector<double> freeEnergyData = thermoProp->getFreeEnergyData();
  const double &fxcr = freeEnergyData[1];
  const double &fxcrr = freeEnergyData[2];
  const double &fxct = freeEnergyData[3];
  const double &fxctt = freeEnergyData[4];
  const double &fxcrt = freeEnergyData[5];
  // QAdder
  const vector<double> QData = thermoProp->getQData();
  const double &Q = QData[0];
  const double &Qr = QData[1];
  const double &Qt = QData[2];
  // Alpha
  double numer = Q - (1.0 / 6.0) * fxcrr + (1.0 / 3.0) * fxcr
                 - (2.0 / 3.0) * (fxctt + fxcrt) + (1.0 / 3.0) * fxct;
  double denom = Q + (1.0 / 3.0) * Qr + (2.0 / 3.0) * Qt;
  return numer / denom;
}

void QVSStls::updateSolution() {
  // Update the structural properties used for output
  ssf = thermoProp->getSsf();
  lfc = thermoProp->getLfc();
}

void QVSStls::init() { Rpa::init(); }

// -----------------------------------------------------------------
// QThermoProp class
// -----------------------------------------------------------------

QThermoProp::QThermoProp(const std::shared_ptr<const QVSStlsInput> &in_)
    : ThermoPropBase(in_),
      structProp(make_shared<QStructProp>(in_)) {
  if (isZeroDegeneracy) {
    throwError(
        "Ground state calculations are not implemented for this scheme.");
  }
  ThermoPropBase::structProp = structProp;
}

vector<double> QThermoProp::getQData() const {
  // QAdder
  const std::vector<double> qVec = structProp->getQ();
  const std::vector<double> rs = structProp->getCouplingParameters();
  const double q = qVec[SIdx::RS_THETA] / rs[SIdx::RS_THETA];
  // QAdder derivative with respect to the coupling parameter
  double qr;
  {
    const double drs = rs[SIdx::RS_UP_THETA] - rs[SIdx::RS_THETA];
    const double &q0 = qVec[SIdx::RS_UP_THETA];
    const double &q1 = qVec[SIdx::RS_DOWN_THETA];
    qr = (q0 - q1) / (2.0 * drs) - q;
  }
  // QAdder derivative with respect to the degeneracy parameter
  double qt;
  {
    const std::vector<double> theta = structProp->getDegeneracyParameters();
    const double dt = theta[SIdx::RS_THETA_UP] - theta[SIdx::RS_THETA];
    const double q0 = qVec[SIdx::RS_THETA_UP] / rs[SIdx::RS_THETA];
    const double q1 = qVec[SIdx::RS_THETA_DOWN] / rs[SIdx::RS_THETA];
    qt = theta[SIdx::RS_THETA] * (q0 - q1) / (2.0 * dt);
  }
  return vector<double>({q, qr, qt});
}

// -----------------------------------------------------------------
// QStructProp class
// -----------------------------------------------------------------

QStructProp::QStructProp(const std::shared_ptr<const QVSStlsInput> &in_)
    : StructPropBase(in_) {
  setupCSR();
  setupCSRDependencies();
}

void QStructProp::setupCSR() {
  std::vector<QVSStlsInput> inVector = setupCSRInput();
  for (const auto &inTmp : inVector) {
    csr.push_back(
        make_shared<QstlsCSR>(make_shared<const QVSStlsInput>(inTmp)));
  }
}

std::vector<QVSStlsInput> QStructProp::setupCSRInput() {
  const double &drs = in().getCouplingResolution();
  const double &dTheta = in().getDegeneracyResolution();
  // If there is a risk of having negative state parameters, shift the
  // parameters so that rs - drs = 0 and/or theta - dtheta = 0
  const double rs = std::max(in().getCoupling(), drs);
  const double theta = std::max(in().getDegeneracy(), dTheta);
  // Setup objects
  std::vector<QVSStlsInput> out;
  for (const double &thetaTmp : {theta - dTheta, theta, theta + dTheta}) {
    for (const double &rsTmp : {rs - drs, rs, rs + drs}) {
      QVSStlsInput inTmp = in();
      inTmp.setDegeneracy(thetaTmp);
      inTmp.setCoupling(rsTmp);
      out.push_back(inTmp);
    }
  }
  // Avoid recomputing the fixed component for the state points with perturbed
  // coupling parameter
  if (in().getFixedRunId() == DEFAULT_INT) {
    for (const int idx : {RS_THETA_DOWN,
                          RS_UP_THETA_DOWN,
                          RS_THETA,
                          RS_UP_THETA,
                          RS_THETA_UP,
                          RS_UP_THETA_UP}) {
      out[idx].setFixedRunId(in().getDatabaseInfo().runId);
    }
  }
  return out;
}

vector<double> QStructProp::getQ() const {
  for (size_t i = 0; i < csr.size(); ++i) {
    const shared_ptr<QstlsCSR> c =
        std::dynamic_pointer_cast<QstlsCSR, CSR>(csr[i]);
    outVector[i] = c->getQAdder();
  }
  return outVector;
}

// -----------------------------------------------------------------
// QstlsCSR class
// -----------------------------------------------------------------

QstlsCSR::QstlsCSR(const std::shared_ptr<const QVSStlsInput> &in_)
    : CSR(),
      Qstls(in_, false),
      itg2D(std::make_shared<Integrator2D>(ItgType::DEFAULT,
                                           in_->getIntError())) {
  if (inRpa().getDegeneracy() == 0.0) {
    throwError("Ground state calculations are not available "
               "for the quantum VS scheme");
  }
  const bool segregatedItg = inRpa().getInt2DScheme() == "segregated";
  if (segregatedItg) { itgGrid = wvg; }
}

void QstlsCSR::init() {
  const string &theory = inRpa().getTheory();
  switch (lfcTheta.type) {
  case CENTERED:
    adrFixedDatabaseName = formatUtil::format("{}_THETA", theory);
    break;
  case FORWARD:
    adrFixedDatabaseName = formatUtil::format("{}_THETA_DOWN", theory);
    break;
  case BACKWARD:
    adrFixedDatabaseName = formatUtil::format("{}_THETA_UP", theory);
    break;
  }
  Qstls::init();
}

void QstlsCSR::computeLfcStls() {
  Qstls::computeLfc();
  *CSR::lfc = Qstls::lfc;
}

void QstlsCSR::computeLfc() {
  Vector2D lfcDerivative = getDerivativeContribution();
  Qstls::lfc.diff(lfcDerivative);
}

double QstlsCSR::getQAdder() const {
  const shared_ptr<Interpolator1D> ssfItp =
      make_shared<Interpolator1D>(wvg, ssf);
  QAdder QTmp(inRpa().getDegeneracy(),
              mu,
              wvg.front(),
              wvg.back(),
              itgGrid,
              itg,
              itg2D,
              ssfItp);
  return QTmp.get();
}

// -----------------------------------------------------------------
// QAdder class
// -----------------------------------------------------------------

// SSF interpolation
double QAdder::ssf(const double &y) const { return interp->eval(y); }

// Denominator integrand
double QAdder::integrandDenominator(const double y) const {
  const double y2 = y * y;
  return 1.0 / (exp(y2 / Theta - mu) + 1.0);
}

// Numerator integrand1
double QAdder::integrandNumerator1(const double q) const {
  const double w = itg2->getX();
  if (q == 0.0) { return 0.0; };
  double w2 = w * w;
  double q2 = q * q;
  double logarg = (w + 2 * q) / (w - 2 * q);
  logarg = (logarg < 0.0) ? -logarg : logarg;
  if (w == 0.0) { return 1.0 / (12.0 * (exp(q2 / Theta - mu) + 1.0)); };
  return q2 / (exp(q2 / Theta - mu) + 1.0) * (q / w * log(logarg) - 1.0) / w2;
}

// Numerator integrand2
double QAdder::integrandNumerator2(const double w) const {
  return (ssf(w) - 1.0);
}

// Denominator integral
void QAdder::getIntDenominator(double &res) const {
  auto func = [&](double y) -> double { return integrandDenominator(y); };
  itg1->compute(func, ItgParam(limits.first, limits.second));
  res = itg1->getSolution();
}

// Get total QAdder
double QAdder::get() const {
  double Denominator;
  getIntDenominator(Denominator);
  auto func1 = [&](const double &w) -> double {
    return integrandNumerator2(w);
  };
  auto func2 = [&](const double &q) -> double {
    return integrandNumerator1(q);
  };
  itg2->compute(
      func1,
      func2,
      Itg2DParam(limits.first, limits.second, limits.first, limits.second),
      itgGrid);
  return 12.0 / (M_PI * numUtil::lambda) * itg2->getSolution() / Denominator;
}
