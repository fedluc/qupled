#include "chemical_potential.hpp"
#include "numerics.hpp"
#include "util.hpp"
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std;
using namespace parallelUtil;

void ChemicalPotential::compute(const vector<double> &guess) {
  auto func = [&](const double &mu) -> double {
    return normalizationCondition(mu);
  };
  BrentRootSolver rsol;
  rsol.solve(func, guess);
  mu = rsol.getSolution();
}

double ChemicalPotential::normalizationCondition(const double &mu) const {
  return gsl_sf_gamma(1.5) * gsl_sf_fermi_dirac_half(mu) -
         2.0 / (3.0 * pow(Theta, 3.0 / 2.0));
}
