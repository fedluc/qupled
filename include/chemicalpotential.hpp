#ifndef CHEMICALPOTENTIAL_HPP
#define CHEMICALPOTENTIAL_HPP

#include "numerics.hpp"

using namespace std;
// using namespace inpututil;

class ChemicalPotential {

private:

  // Degeneracy parameter
  const double Theta;
  // Chemical potential
  double mu;
  // Normalization condition
  double normalizationCondition(double mu) const;
  
public:

  ChemicalPotential(double Theta_) : Theta(Theta_) {};
  void compute(const vector<double> &guess);
  double get() const {return mu;};
  
};

#endif
