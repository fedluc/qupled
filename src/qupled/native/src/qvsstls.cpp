#include "qvsstls.hpp"
#include "format.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;
using namespace MPIUtil;
using namespace GridPoints;
using ItgType   = Integrator1D::Type;
using ItgParam  = Integrator1D::Param;
using Itg2DParam = Integrator2D::Param;

// -----------------------------------------------------------------
// VSQstlsWorker
// -----------------------------------------------------------------

VSQstlsWorker::VSQstlsWorker(const std::shared_ptr<const QVSStlsInput> &in,
                             bool verbose,
                             GridPoint p)
    : Qstls(in, verbose) {
  // Set adrFixedDatabaseName based on theta offset (mirrors QstlsCSR::initWorker)
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

double VSQstlsWorker::computeQAdder(
    const std::shared_ptr<Integrator2D> &itg2D,
    const std::vector<double> &itgGrid) const {
  const auto ssfItp = make_shared<Interpolator1D>(wvg, ssf);
  QAdder q(inPtr->getDegeneracy(),
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
// QVSStls
// -----------------------------------------------------------------

QVSStls::QVSStls(const std::shared_ptr<const QVSStlsInput> &in)
    : VSBase(),
      inPtr(in),
      grid(in),
      itg2D(make_shared<Integrator2D>(ItgType::DEFAULT, in->getIntError())) {
  if (in->getDegeneracy() == 0.0) {
    throwError(
        "Ground state calculations are not implemented for this scheme.");
  }
  const bool segregatedItg = in->getInt2DScheme() == "segregated";
  if (segregatedItg) { itgGrid = grid.centralWorker().getWvg(); }
  setRsGrid();
  setFxcIntegrand();
}

const VSInput &QVSStls::in() const {
  return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
}

const Input &QVSStls::inScheme() const { return *inPtr; }

void QVSStls::init() {
  // Worker initialisation is deferred to StatePointGrid::compute()
}

void QVSStls::updateSolution() {
  const GridPoint out = getOutputGridPoint();
  ssf = grid.getSsf(out);
  lfc = grid.getLfc(out);
}

int QVSStls::runGrid() {
  grid.setAlpha(alpha);
  return grid.compute();
}

double QVSStls::getCoupling(GridPoint p) const { return grid.getCoupling(p); }

double QVSStls::getDegeneracy(GridPoint p) const {
  return grid.getDegeneracy(p);
}

double QVSStls::getFxcIntegrandValue(GridPoint p) const {
  return grid.getFxcIntegrandValue(p);
}

vector<double> QVSStls::computeQData() {
  const double rs = grid.getCoupling(CENTER);
  // Q at the target state point
  const double q =
      grid.centralWorker().computeQAdder(itg2D, itgGrid) / rs;
  // Derivative w.r.t. coupling
  double qr;
  {
    const double drs = grid.getCoupling(RS_UP_THETA) - rs;
    const auto getQ  = [&](GridPoint p) {
      return grid.getWorkerAt(p).computeQAdder(itg2D, itgGrid);
    };
    const double q0 = getQ(RS_UP_THETA);
    const double q1 = getQ(RS_DOWN_THETA);
    qr = (q0 - q1) / (2.0 * drs) - q;
  }
  // Derivative w.r.t. degeneracy
  double qt;
  {
    const double theta = grid.getDegeneracy(CENTER);
    const double dt    = grid.getDegeneracy(RS_THETA_UP) - theta;
    const double q0 = grid.getWorkerAt(RS_THETA_UP).computeQAdder(itg2D, itgGrid) / rs;
    const double q1 = grid.getWorkerAt(RS_THETA_DOWN).computeQAdder(itg2D, itgGrid) / rs;
    qt = theta * (q0 - q1) / (2.0 * dt);
  }
  return {q, qr, qt};
}

// Delegation to central worker for Python-exposed properties

const vector<double> &QVSStls::getWvg() const {
  return grid.centralWorker().getWvg();
}

const Vector2D &QVSStls::getIdr() const {
  return grid.centralWorker().getIdr();
}

vector<double> QVSStls::getSdr() const {
  return grid.centralWorker().getSdr();
}

double QVSStls::getUInt() const { return grid.centralWorker().getUInt(); }

double QVSStls::getError() const { return grid.getError(); }

// -----------------------------------------------------------------
// QAdder
// -----------------------------------------------------------------

double QAdder::ssf(const double &y) const { return interp->eval(y); }

double QAdder::integrandDenominator(const double y) const {
  const double y2 = y * y;
  return 1.0 / (exp(y2 / Theta - mu) + 1.0);
}

double QAdder::integrandNumerator1(const double q) const {
  const double w = itg2->getX();
  if (q == 0.0) { return 0.0; }
  const double w2  = w * w;
  const double q2  = q * q;
  double logarg = (w + 2 * q) / (w - 2 * q);
  logarg = (logarg < 0.0) ? -logarg : logarg;
  if (w == 0.0) { return 1.0 / (12.0 * (exp(q2 / Theta - mu) + 1.0)); }
  return q2 / (exp(q2 / Theta - mu) + 1.0) * (q / w * log(logarg) - 1.0) / w2;
}

double QAdder::integrandNumerator2(const double w) const {
  return (ssf(w) - 1.0);
}

void QAdder::getIntDenominator(double &res) const {
  auto func = [&](double y) -> double { return integrandDenominator(y); };
  itg1->compute(func, ItgParam(limits.first, limits.second));
  res = itg1->getSolution();
}

double QAdder::get() const {
  double Denominator;
  getIntDenominator(Denominator);
  auto func1 = [&](const double &w) -> double { return integrandNumerator2(w); };
  auto func2 = [&](const double &q) -> double { return integrandNumerator1(q); };
  itg2->compute(
      func1,
      func2,
      Itg2DParam(limits.first, limits.second, limits.first, limits.second),
      itgGrid);
  return 12.0 / (M_PI * numUtil::lambda) * itg2->getSolution() / Denominator;
}
