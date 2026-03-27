#include "vs/qvsstls.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;
using namespace MPIUtil;
using namespace GridPoints;
using ItgType = Integrator1D::Type;
using ItgParam = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;

// -----------------------------------------------------------------
// VSQstlsWorker
// -----------------------------------------------------------------

VSQstlsWorker::VSQstlsWorker(const std::shared_ptr<const QVSStlsInput> &in,
                             GridPoint p)
    : Qstls(in, false) {
  const string &theory = in->getTheory();
  switch (p.theta) {
  case GridPoint::Theta::DOWN:
    adrFixedDatabaseName = formatUtil::format("{}_THETA_DOWN", theory);
    break;
  case GridPoint::Theta::CENTER:
    adrFixedDatabaseName = formatUtil::format("{}_THETA", theory);
    break;
  case GridPoint::Theta::UP:
    adrFixedDatabaseName = formatUtil::format("{}_THETA_UP", theory);
    break;
  }
}

double VSQstlsWorker::getQAdder() const {
  const bool segregatedItg = inPtr->getInt2DScheme() == "segregated";
  const vector<double> itgGrid = (segregatedItg) ? wvg : vector<double>();
  shared_ptr<Integrator2D> itg2D =
      make_shared<Integrator2D>(inPtr->getIntError());
  const auto ssfItp = make_shared<Interpolator1D>(wvg, ssf);
  const QAdder q = QAdder(inPtr->getDegeneracy(),
                          mu,
                          wvg.front(),
                          wvg.back(),
                          itgGrid,
                          itg,
                          itg2D,
                          ssfItp);
  return q.get();
}

// -----------------------------------------------------------------
// VSQstlsManager
// -----------------------------------------------------------------

VSQstlsManager::VSQstlsManager(const std::shared_ptr<const QVSStlsInput> &in)
    : VSManager(),
      Qstls(in, false),
      managerInPtr_(in),
      itg2D(make_shared<Integrator2D>(ItgType::DEFAULT, in->getIntError())) {
  const double drs_ = in->getCouplingResolution();
  const double dTheta_ = in->getDegeneracyResolution();
  const double rs0 = std::max(in->getCoupling(), drs_);
  const double theta0 = std::max(in->getDegeneracy(), dTheta_);
  for (const auto tOff : {GridPoint::Theta::DOWN,
                          GridPoint::Theta::CENTER,
                          GridPoint::Theta::UP}) {
    const double thetaTmp = theta0 + static_cast<int>(tOff) * dTheta_;
    for (const auto rOff :
         {GridPoint::Rs::DOWN, GridPoint::Rs::CENTER, GridPoint::Rs::UP}) {
      const double rsTmp = rs0 + static_cast<int>(rOff) * drs_;
      auto inTmp = std::make_shared<QVSStlsInput>(*in);
      inTmp->setCoupling(rsTmp);
      inTmp->setDegeneracy(thetaTmp);
      if (rOff != GridPoint::Rs::DOWN && in->getFixedRunId() == DEFAULT_INT) {
        inTmp->setFixedRunId(in->getFixedRunId());
      }
      const GridPoint gp{rOff, tOff};
      const size_t idx = gp.toIndex();
      rsValues[idx] = rsTmp;
      thetaValues[idx] = thetaTmp;
      workers[idx] = std::make_unique<VSQstlsWorker>(inTmp, gp);
    }
  }
  setupDerivativeData();
  // Setup integration grid
  const bool segregatedItg = in->getInt2DScheme() == "segregated";
  if (segregatedItg) { itgGrid = VSManager::getWvg(GridPoints::CENTER); }
}

// -----------------------------------------------------------------
// QVSStls
// -----------------------------------------------------------------

QVSStls::QVSStls(const std::shared_ptr<const QVSStlsInput> &in)
    : VSBase(),
      inPtr(in),
      grid_(in) {
  if (in->getDegeneracy() == 0.0) {
    throwError(
        "Ground state calculations are not implemented for this scheme.");
  }
  setRsGrid();
  setFxcIntegrand();
}

const VSInput &QVSStls::in() const {
  return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
}

const Input &QVSStls::inScheme() const { return *inPtr; }

// -----------------------------------------------------------------
// QAdder
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