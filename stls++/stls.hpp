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
  // Compute Hartree-Fock static structure factor
  void computeSsfHF();

public:

  // Constructor
  Stls(Input in_) : in(in_), computedChemicalPotential(false) {;};
  // Compute stls scheme
  void compute(); 


};

// Class for the ideal density response calculation
class Idr {

private:
  
  // Matsubara frequency
  const int l = 0;
  // Frequency for zero temperature calculations
  const double Omega = 0;
  // Wave-vector
  const double x = 0;
  // Degeneracy parameter
  const double Theta = 0;
  // Chemical potential
  const double mu = 0;
  // Integration limits for finite temperature calculations
  const double yMin = 0;
  const double yMax = 0;
  // Idr integrand for frequency = l and wave-vector x
  double xl(double y);
  // Idr integrand for frequency = 0 and wave-vector x
  double x0(double y);  
  
public:

  // Constructor for finite temperature calculations
  Idr(int l_, double x_, double Theta_,
      double mu_, double yMin_, double yMax_)
    : l(l_), x(x_), Theta(Theta_), mu(mu_),
      yMin(yMin_), yMax(yMax_) {;};
  // Constructor for zero temperature calculations
  Idr(double Omega_, double x_) : Omega(Omega_), x(x_) {;};
  // Get at finite temperature
  double get(Integrator1D &itg);
  // Integrand
  double integrand(double y);
  // Get real part at zero temperature
  double re0();
  // Get imaginary part at zero temperature
  double im0();
  // Get frequency derivative of the real part at zero temperature
  double re0Der();
  
};



// Class for the Hartree-Fock static structure factor
class SsfHF {

private:
  
  // Wave-vector
  const double x = 0;
  // Degeneracy parameter
  const double Theta = 0;
  // Chemical potential
  const double mu = 0;
  // Integration limits for finite temperature calculations
  const double yMin = 0;
  const double yMax = 0;
  
public:

  // Constructor for finite temperature calculations
  SsfHF(double x_, double Theta_, double mu_, double yMin_, double yMax_)
    : x(x_), Theta(Theta_), mu(mu_), yMin(yMin_), yMax(yMax_) {;};
  // Constructor for zero temperature calculations
  SsfHF(double x_) : x(x_) {;};
  // Get at finite temperature
  double get(Integrator1D &itg);
  // Get integrand
  double integrand(double y);
  // Get at zero temperature
  double get0();
  
};
 
#endif
