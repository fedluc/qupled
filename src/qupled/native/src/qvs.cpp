#include "qvs.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"
#include <filesystem>
#include <fmt/core.h>

using namespace std;
using namespace vecUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;
using ItgType = Integrator1D::Type;

// -----------------------------------------------------------------
// QVSStls class
// -----------------------------------------------------------------

QVSStls::QVSStls(const QVSStlsInput &in_)
    : VSBase(in_),
      Qstls(in_, false, false),
      in(in_),
      thermoProp(make_shared<QThermoProp>(in_)) {
  VSBase::thermoProp = thermoProp;
}

QVSStls::QVSStls(const QVSStlsInput &in_, const QThermoProp &thermoProp_)
    : VSBase(in_, false),
      Qstls(in_, false, false),
      in(in_),
      thermoProp(make_shared<QThermoProp>(in_)) {
  VSBase::thermoProp = thermoProp;
  thermoProp->copyFreeEnergyIntegrand(thermoProp_);
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
  double numer = Q - (1.0 / 6.0) * fxcrr + (1.0 / 3.0) * fxcr;
  double denom = Q + (1.0 / 3.0) * Qr;
  if (in.getDegeneracy() > 0.0) {
    numer += -(2.0 / 3.0) * fxctt - (2.0 / 3.0) * fxcrt + (1.0 / 3.0) * fxct;
    denom += (2.0 / 3.0) * Qt;
  }
  return numer / denom;
}

void QVSStls::updateSolution() {
  // Update the structural properties used for output
  adr = thermoProp->getAdr();
  ssf = thermoProp->getSsf();
  slfc = thermoProp->getSlfc();
}

void QVSStls::initScheme() { Rpa::init(); }

void QVSStls::initFreeEnergyIntegrand() {
  if (!thermoProp->isFreeEnergyIntegrandIncomplete()) { return; }
  println("Missing points in the free energy integrand: subcalls will be "
          "performed to collect the necessary data");
  println("-----------------------------------------------------------------"
          "----------");
  QVSStlsInput inTmp = in;
  while (thermoProp->isFreeEnergyIntegrandIncomplete()) {
    const double rs = thermoProp->getFirstUnsolvedStatePoint();
    println(fmt::format("Subcall: solving qVS scheme for rs = {:.5f}", rs));
    inTmp.setCoupling(rs);
    QVSStls scheme(inTmp, *thermoProp);
    scheme.compute();
    thermoProp->copyFreeEnergyIntegrand(*(scheme.thermoProp));
    println("Done");
    println("-----------------------------------------------------------------"
            "----------");
  }
  println("Subcalls completed");
}

// -----------------------------------------------------------------
// QThermoProp class
// -----------------------------------------------------------------

QThermoProp::QThermoProp(const QVSStlsInput &in_)
    : ThermoPropBase(in_, in_),
      structProp(make_shared<QStructProp>(in_)) {
  if (isZeroDegeneracy) {
    throwError("Ground state calculations are not available "
               "for the quantum VS scheme");
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

const Vector2D &QThermoProp::getAdr() {
  if (!structProp->isComputed()) { structProp->compute(); }
  return structProp->getCsr(getStructPropIdx()).getAdr();
}

// -----------------------------------------------------------------
// QStructProp class
// -----------------------------------------------------------------

QStructProp::QStructProp(const QVSStlsInput &in_)
    : Logger(MPIUtil::isRoot()),
      StructPropBase(),
      in(in_) {
  setupCSR();
  setupCSRDependencies();
}

void QStructProp::setupCSR() {
  std::vector<QVSStlsInput> inVector = setupCSRInput();
  for (const auto &inTmp : inVector) {
    csr.push_back(make_shared<QstlsCSR>(inTmp));
  }
  for (const auto &c : csr) {
    StructPropBase::csr.push_back(c);
  }
}

std::vector<QVSStlsInput> QStructProp::setupCSRInput() {
  const double &drs = in.getCouplingResolution();
  const double &dTheta = in.getDegeneracyResolution();
  // If there is a risk of having negative state parameters, shift the
  // parameters so that rs - drs = 0 and/or theta - dtheta = 0
  const double rs = std::max(in.getCoupling(), drs);
  const double theta = std::max(in.getDegeneracy(), dTheta);
  // Setup objects
  std::vector<QVSStlsInput> out;
  for (const double &thetaTmp : {theta - dTheta, theta, theta + dTheta}) {
    for (const double &rsTmp : {rs - drs, rs, rs + drs}) {
      QVSStlsInput inTmp = in;
      inTmp.setDegeneracy(thetaTmp);
      inTmp.setCoupling(rsTmp);
      out.push_back(inTmp);
    }
  }
  return out;
}

const QstlsCSR &QStructProp::getCsr(const Idx &idx) const { return *csr[idx]; }

void QStructProp::doIterations() {
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
#pragma omp parallel num_threads(ompThreads) if (useOMP)
    {
#pragma omp for
      for (auto &c : csr) {
        c->computeAdrQStls();
      }
#pragma omp for
      for (size_t i = 0; i < csr.size(); ++i) {
        auto &c = csr[i];
        c->computeAdr();
        c->computeSsf();
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
  // Set static structure factor for output
  for (auto &c : csr) {
    c->updateSsf();
  }
}

vector<double> QStructProp::getQ() const {
  for (size_t i = 0; i < csr.size(); ++i) {
    outVector[i] = csr[i]->getQAdder();
  }
  return outVector;
}

// -----------------------------------------------------------------
// QstlsCSR class
// -----------------------------------------------------------------

QstlsCSR::QstlsCSR(const QVSStlsInput &in_)
    : CSR(in_, in_),
      Qstls(in_, false, false),
      in(in_) {
  if (in.getDegeneracy() == 0.0) {
    throwError("Ground state calculations are not available "
               "for the quantum VS scheme");
  }
}

void QstlsCSR::init() {
  switch (lfcTheta.type) {
  case CENTERED: adrFixedFileName = "THETA.bin"; break;
  case FORWARD: adrFixedFileName = "THETA_DOWN.bin"; break;
  case BACKWARD: adrFixedFileName = "THETA_UP.bin"; break;
  }
  if (!in.getFixed().empty()) {
    std::filesystem::path fullPath = in.getFixed();
    fullPath /= adrFixedFileName;
    adrFixedFileName = fullPath.string();
  }
  if (std::filesystem::exists(adrFixedFileName)) {
    Stls::init();
    readAdrFixedFile(adrFixed, adrFixedFileName, false);
  } else {
    Qstls::init();
  }
  // MPI barrier to make sure that all processes see the same files
  MPIUtil::barrier();
}

void QstlsCSR::computeAdrQStls() {
  Qstls::computeAdr();
  *lfc = adr;
}

Vector2D QstlsCSR::getDerivativeContribution() const {
  Vector2D out = CSR::getDerivativeContribution();
  out.linearCombination(*lfc, -alpha / 3.0);
  return out;
}

void QstlsCSR::computeAdr() {
  Vector2D adrDerivative = getDerivativeContribution();
  adr.diff(adrDerivative);
}

double QstlsCSR::getQAdder() const {
  Integrator1D itg1(ItgType::DEFAULT, in.getIntError());
  Integrator2D itg2(ItgType::DEFAULT, ItgType::DEFAULT, in.getIntError());
  const bool segregatedItg = in.getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  const Interpolator1D ssfItp(wvg, ssf);
  QAdder QTmp(in.getDegeneracy(),
              mu,
              wvg.front(),
              wvg.back(),
              itgGrid,
              itg1,
              itg2,
              ssfItp);
  return QTmp.get();
}

// -----------------------------------------------------------------
// QAdder class
// -----------------------------------------------------------------

// SSF interpolation
double QAdder::ssf(const double &y) const { return interp.eval(y); }

// Denominator integrand
double QAdder::integrandDenominator(const double y) const {
  const double y2 = y * y;
  return 1.0 / (exp(y2 / Theta - mu) + 1.0);
}

// Numerator integrand1
double QAdder::integrandNumerator1(const double q) const {
  const double w = itg2.getX();
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
  itg1.compute(func, ItgParam(limits.first, limits.second));
  res = itg1.getSolution();
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
  itg2.compute(
      func1,
      func2,
      Itg2DParam(limits.first, limits.second, limits.first, limits.second),
      itgGrid);
  return 12.0 / (M_PI * lambda) * itg2.getSolution() / Denominator;
}
