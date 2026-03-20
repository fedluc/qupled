#include "vsstls.hpp"
#include "input.hpp"
#include "mpi_util.hpp"
#include "thermo_util.hpp"

using namespace std;
using namespace MPIUtil;
using namespace GridPoints;

// -----------------------------------------------------------------
// VSStls
// -----------------------------------------------------------------

VSStls::VSStls(const std::shared_ptr<const VSStlsInput> &in)
    : VSBase(),
      inPtr(in),
      grid(in) {
  if (in->getDimension() == dimensionsUtil::Dimension::D2) {
    throwError("2D calculations are not implemented for this scheme.");
  }
  setRsGrid();
  setFxcIntegrand();
}

const VSInput &VSStls::in() const {
  return *StlsUtil::dynamic_pointer_cast<Input, VSInput>(inPtr);
}

const Input &VSStls::inScheme() const { return *inPtr; }

void VSStls::init() {
  // Worker initialisation is deferred to StatePointGrid::compute()
}

void VSStls::updateSolution() {
  const GridPoint out = getOutputGridPoint();
  ssf = grid.getSsf(out);
  lfc = grid.getLfc(out);
}

int VSStls::runGrid() {
  grid.setAlpha(alpha);
  return grid.compute();
}

double VSStls::getCoupling(GridPoint p) const { return grid.getCoupling(p); }

double VSStls::getDegeneracy(GridPoint p) const {
  return grid.getDegeneracy(p);
}

double VSStls::getFxcIntegrandValue(GridPoint p) const {
  return grid.getFxcIntegrandValue(p);
}

vector<double> VSStls::computeQData() {
  const double u = grid.getUInt(CENTER);
  // Derivative w.r.t. coupling
  double ur;
  {
    const double drs = grid.getCoupling(RS_UP_THETA)
                       - grid.getCoupling(CENTER);
    const double u0 = grid.getFxcIntegrandValue(RS_UP_THETA);
    const double u1 = grid.getFxcIntegrandValue(RS_DOWN_THETA);
    ur = (u0 - u1) / (2.0 * drs) - u;
  }
  // Derivative w.r.t. degeneracy
  double ut;
  {
    const double theta = grid.getDegeneracy(CENTER);
    const double dt    = grid.getDegeneracy(RS_THETA_UP) - theta;
    const double u0    = grid.getUInt(RS_THETA_UP);
    const double u1    = grid.getUInt(RS_THETA_DOWN);
    ut = theta * (u0 - u1) / (2.0 * dt);
  }
  return {u, ur, ut};
}

// Delegation to central worker for Python-exposed properties

const vector<double> &VSStls::getWvg() const {
  return grid.centralWorker().getWvg();
}

const Vector2D &VSStls::getIdr() const {
  return grid.centralWorker().getIdr();
}

vector<double> VSStls::getSdr() const {
  return grid.centralWorker().getSdr();
}

double VSStls::getUInt() const { return grid.centralWorker().getUInt(); }

double VSStls::getError() const { return grid.getError(); }
