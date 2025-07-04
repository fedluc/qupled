#include "chemical_potential.hpp"
#include "numerics.hpp"

using namespace std;

void ChemicalPotential::compute(const vector<double> &guess) {
  auto func = [&](const double &mu) -> double {
    return normalizationCondition(mu);
  };
  BrentRootSolver rsol;
  rsol.solve(func, guess);
  mu = rsol.getSolution();
}

double ChemicalPotential::normalizationCondition(const double &mu) const {
  return SpecialFunctions::fermiDirac(mu) - 2.0 / (3.0 * pow(Theta, 1.5));
}

void ChemicalPotential::compute2D() {
  mu = log(exp(1/Theta) - 1.0);
}