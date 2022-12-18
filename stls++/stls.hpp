#ifndef STLS_HPP
#define STLS_HPP

#include <vector>
#include "input.hpp"
#include "numerics.hpp"

using namespace std;

class Stls {

private: 

  // Wave vector grid
  vector<double> wvg;
  // Ideal density response
  vector<vector<double>> idr;
  // Static local field correction
  vector<double> slfc;
  // Static structure factor
  vector<double> ssf;
  // Hartree-Fock static structure factor
  vector<double> ssfHF;
  // Input data
  Input in;
  // Chemical potential
  double mu;
  bool computedChemicalPotential;
  // Construct wave vector grid
  void buildWvGrid();
  // Compute chemical potential
  void computeChemicalPotential();
  // Compute the ideal density response
  void computeIdr();
  void computeIdrSingleFrequency(Integrator1D &itg,
				 const int l);

public:

  // Constructor
  Stls(Input in_) : in(in_), computedChemicalPotential(false) {;};
  // Compute stls scheme
  void compute(); 


};

// Class for the ideal density response integration
class Idr {

private:
  
  // Matsubara frequency
  const int l;
  // Wave-vector
  const double x;
  // Degeneracy parameter
  const double Theta;
  // Chemical potential
  const double mu;
  // Idr integrand for frequency = l and wave-vector x
  double xl(double y);
  // Idr integrand for frequency = 0 and wave-vector x
  double x0(double y);
  
public:

  // Constructor
  Idr(int l_, double x_, double Theta_, double mu_)
    : l(l_), x(x_), Theta(Theta_), mu(mu_) {;};
  // Get integrand
  double get(double y);

};

#endif
