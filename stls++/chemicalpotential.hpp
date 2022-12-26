#ifndef CHEMICALPOTENTIAL_HPP
#define CHEMICALPOTENTIAL_HPP

#include <vector>
#include <iostream>
#include "inpututil.hpp"

using namespace std;
using namespace inpututil;

class ChemicalPotential {

private:

  double normalizationCondition(double mu);
  double Theta;
  double mu;
  
public:

  ChemicalPotential(double Theta_) : Theta(Theta_) {};
  void compute(cVector<double> &guess);
  double get() {return mu;};
  
};

#endif
