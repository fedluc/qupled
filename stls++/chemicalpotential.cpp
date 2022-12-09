#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "chemicalpotential.hpp"
#include "numerics.hpp"

void ChemicalPotential::compute(const vector<double> &guess){
  auto func = [this](double mu)->double{return normalizationCondition(mu);};
  RootSolver rs;
  rs.solve(func, guess);
  mu = rs.getSolution();
}

double  ChemicalPotential::normalizationCondition(double mu){
  return gsl_sf_gamma(1.5)*gsl_sf_fermi_dirac_half(mu) 
         - 2.0/(3.0*pow(Theta, 3.0/2.0));
}
