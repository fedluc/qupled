#include "chemical_potential.hpp"
#include "dimensions_util.hpp"
#include "numerics.hpp"

using namespace std;

void ChemicalPotential::compute3D() {
  auto func = [&](const double &mu) -> double {
    return normalizationCondition(mu);
  };
  BrentRootSolver rsol;
  rsol.solve(func, in->getChemicalPotentialGuess());
  mu = rsol.getSolution();
}

void ChemicalPotential::compute2D() {
  const double &Theta = in->getDegeneracy();
  mu = log(exp(1.0 / Theta) - 1.0);
}

double ChemicalPotential::normalizationCondition(const double &mu) const {
  const double &Theta = in->getDegeneracy();
  return SpecialFunctions::fermiDirac(mu) - 2.0 / (3.0 * pow(Theta, 1.5));
}