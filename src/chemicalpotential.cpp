#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "chemicalpotential.hpp"

void ChemicalPotential::compute(const vector<double> &guess){
  auto func = [this](double mu)->double{return normalizationCondition(mu);};
  RootSolver rsol;
  rsol.solve(func, guess);
  if (!rsol.success()) {
    throw runtime_error("Chemical potential: the root solver "
			"did not converge to the desired accuracy.");
  }
  mu = rsol.getSolution();
}

double  ChemicalPotential::normalizationCondition(double mu) const {
  return gsl_sf_gamma(1.5)*gsl_sf_fermi_dirac_half(mu) 
         - 2.0/(3.0*pow(Theta, 3.0/2.0));
}
