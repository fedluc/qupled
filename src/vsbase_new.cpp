#include "vsbase_new.hpp"
#include "input.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector_util.hpp"

using namespace std;

// -----------------------------------------------------------------
// VSBase class
// -----------------------------------------------------------------

int VSBase::compute() {
  try {
    init();
    if (verbose) cout << "Free parameter calculation ..." << endl;
    doIterations();
    if (verbose) cout << "Done" << endl;
    return 0;
  } catch (const runtime_error &err) {
    cerr << err.what() << endl;
    return 1;
  }
}

vector<vector<double>> VSBase::getFreeEnergyIntegrand() const {
  return thermoProp->getFreeEnergyIntegrand();
}

vector<double> VSBase::getFreeEnergyGrid() const {
  return thermoProp->getFreeEnergyGrid();
}

vector<double> VSBase::getAlpha() const {
  return thermoProp->getAlpha();
}

void VSBase::doIterations() {
  auto func = [this](const double &alphaTmp) -> double {
    return alphaDifference(alphaTmp);
  };
  SecantSolver rsol(in.getErrMinAlpha(), in.getNIterAlpha());
  rsol.solve(func, in.getAlphaGuess());
  alpha = rsol.getSolution();
  if (verbose) { std::cout << "Free parameter = " << alpha << std::endl; }
  updateSolution();
}

double VSBase::alphaDifference(const double &alphaTmp) {
  alpha = alphaTmp;
  thermoProp->setAlpha(alpha);
  const double alphaTheoretical = computeAlpha();
  return alpha - alphaTheoretical;
}

// -----------------------------------------------------------------
// ThermoPropBase class
// -----------------------------------------------------------------

// Add methods from ThermoPropBase here instead of having them in the header
// The recursive calls maybe should be moved in vsbase
