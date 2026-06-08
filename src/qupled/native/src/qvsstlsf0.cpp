#include "qvsstlsf0.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "vector_util.hpp"

using namespace std;
using namespace vecUtil;
using namespace MPIUtil;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;
using ItgType = Integrator1D::Type;

// -----------------------------------------------------------------
// QVSStlsF0 class
// -----------------------------------------------------------------

QVSStlsF0::QVSStlsF0(const std::shared_ptr<const QVSStlsF0Input> &in_)
    : VSBase(),
      QstlsF0(in_),
      thermoProp(make_shared<QThermoPropF0>(in_)) {
  VSBase::thermoProp = thermoProp;
}

double QVSStlsF0::computeAlpha() {
  thermoProp->compute();

  const vector<double> freeEnergyData = thermoProp->getFreeEnergyData();
  const double &fxcr = freeEnergyData[1];
  const double &fxcrr = freeEnergyData[2];
  const double &fxct = freeEnergyData[3];
  const double &fxctt = freeEnergyData[4];
  const double &fxcrt = freeEnergyData[5];

  const vector<double> QData = thermoProp->getQData();
  const double &Q = QData[0];
  const double &Qr = QData[1];
  const double &Qt = QData[2];

  const double numer = Q - (1.0 / 6.0) * fxcrr + (1.0 / 3.0) * fxcr
                       - (2.0 / 3.0) * (fxctt + fxcrt) + (1.0 / 3.0) * fxct;
  const double denom = Q + (1.0 / 3.0) * Qr + (2.0 / 3.0) * Qt;
  return numer / denom;
}

void QVSStlsF0::updateSolution() {
  ssf = thermoProp->getSsf();
  lfc = thermoProp->getLfc();
}

void QVSStlsF0::init() { Rpa::init(); }

// -----------------------------------------------------------------
// QThermoPropF0 class
// -----------------------------------------------------------------

QThermoPropF0::QThermoPropF0(const std::shared_ptr<const QVSStlsF0Input> &in_)
    : ThermoPropBase(in_),
      structProp(make_shared<QstlsF0CSR>(in_)) {
  if (isZeroDegeneracy) {
    throwError("Ground state calculations are not implemented for this scheme.");
  }
  ThermoPropBase::structProp = structProp;
}

vector<double> QThermoPropF0::getQData() const {
  const double q = structProp->getQAdder(SIdx::RS_THETA)
                   / structProp->getCoupling(SIdx::RS_THETA);

  double qr;
  {
    const double drs = structProp->getCoupling(SIdx::RS_UP_THETA)
                       - structProp->getCoupling(SIdx::RS_THETA);
    const double &q0 = structProp->getQAdder(SIdx::RS_UP_THETA);
    const double &q1 = structProp->getQAdder(SIdx::RS_DOWN_THETA);
    qr = (q0 - q1) / (2.0 * drs) - q;
  }

  double qt;
  {
    const double &rs = structProp->getCoupling(SIdx::RS_THETA);
    const double &theta = structProp->getDegeneracy(SIdx::RS_THETA);
    const double dt = structProp->getDegeneracy(SIdx::RS_THETA_UP)
                      - structProp->getDegeneracy(SIdx::RS_THETA);
    const double q0 = structProp->getQAdder(SIdx::RS_THETA_UP) / rs;
    const double q1 = structProp->getQAdder(SIdx::RS_THETA_DOWN) / rs;
    qt = theta * (q0 - q1) / (2.0 * dt);
  }
  return vector<double>({q, qr, qt});
}

// -----------------------------------------------------------------
// QstlsF0CSR class
// -----------------------------------------------------------------

QstlsF0CSR::QstlsF0CSR(const std::shared_ptr<const QVSStlsF0Input> &in_,
                       const bool isMaster_)
    : CSR(isMaster_),
      QstlsF0(in_),
      itg2D(std::make_shared<Integrator2D>(ItgType::DEFAULT,
                                           in_->getIntError())) {
  if (inRpa().getDegeneracy() == 0.0) {
    throwError("Ground state calculations are not available "
               "for the quantum VS scheme");
  }
  const bool segregatedItg = inRpa().getInt2DScheme() == "segregated";
  if (segregatedItg) { itgGrid = wvg; }
  if (isManager) { setupWorkers(*in_); }
}

void QstlsF0CSR::initWorker() {
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
  QstlsF0::init();
}

double QstlsF0CSR::getQAdder(const size_t &idx) const {
  if (isManager) {
    const CSR &baseWorker = *workers[idx];
    const QstlsF0CSR &thisWorker = static_cast<const QstlsF0CSR &>(baseWorker);
    return thisWorker.getQAdder(idx);
  }
  const shared_ptr<Interpolator1D> ssfItp =
      make_shared<Interpolator1D>(wvg, ssf);
  QAdderF0 QTmp(inRpa().getCoupling(),
                inRpa().getDegeneracy(),
                mu,
                in().getPimcEta(),
                in().getPimcYSec(),
                in().getPimcACutoff(),
                wvg.front(),
                wvg.back(),
                itgGrid,
                itg,
                itg2D,
                ssfItp);
  return QTmp.get();
}

// -----------------------------------------------------------------
// QAdderF0 class
// -----------------------------------------------------------------

double QAdderF0::ssf(const double &y) const { return interp->eval(y); }

double QAdderF0::integrandDenominator(const double y) const {
  return f0.value(y);
}

double QAdderF0::integrandNumerator1(const double q) const {
  const double w = itg2->getX();
  if (q == 0.0) { return 0.0; }
  const double w2 = w * w;
  const double q2 = q * q;
  double logarg = (w + 2.0 * q) / (w - 2.0 * q);
  logarg = (logarg < 0.0) ? -logarg : logarg;
  if (w == 0.0) { return f0.value(q) / 12.0; }
  return q2 * f0.value(q) * (q / w * log(logarg) - 1.0) / w2;
}

double QAdderF0::integrandNumerator2(const double w) const {
  return (ssf(w) - 1.0);
}

void QAdderF0::getIntDenominator(double &res) const {
  auto func = [&](double y) -> double { return integrandDenominator(y); };
  itg1->compute(func, ItgParam(limits.first, limits.second));
  res = itg1->getSolution();
}

double QAdderF0::get() const {
  double denominator;
  getIntDenominator(denominator);

  auto func1 = [&](const double &w) -> double { return integrandNumerator2(w); };
  auto func2 = [&](const double &q) -> double { return integrandNumerator1(q); };
  itg2->compute(func1,
                func2,
                Itg2DParam(limits.first, limits.second, limits.first, limits.second),
                itgGrid);
  return 12.0 / (M_PI * numUtil::lambda) * itg2->getSolution() / denominator;
}
