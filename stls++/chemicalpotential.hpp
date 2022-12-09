#ifndef CHEMICALPOTENTIAL_HPP
#define CHEMICALPOTENTIAL_HPP

#include <vector>
#include <iostream>
#include <gsl/gsl_roots.h>
#include "numerics.hpp"

using namespace std;

class ChemicalPotential {

private:

  double normalizationCondition(double mu);
  double Theta;
  double mu;
  
public:

  ChemicalPotential(double Theta_) : Theta(Theta_) {};
  void compute(const vector<double> &guess);
  double get() {return mu;};
  
};

#endif
