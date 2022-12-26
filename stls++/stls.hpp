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
  vector<double> slfcOld;
  vector<double> slfc;
  // Static structure factor
  vector<double> ssf;
  // Hartree-Fock static structure factor
  vector<double> ssfHF;
  // Constant for unit conversion
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  // Integrator object
  shared_ptr<Integrator1D> itg;
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
  void computeIdrSingleFrequency(const int l);
  // Compute Hartree-Fock static structure factor
  void computeSsfHF();
  // Compute static structure factor at finite temperature
  void computeSsf();
  // Compute static local field correction
  void computeSlfc();
  // Iterations to solve the stls scheme
  void doIterations();
  void initialGuess();
  double computeError();
  void updateSolution();

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
  // Idr integrand for frequency = l and wave-vector y
  double xl(const double y);
  // Idr integrand for frequency = 0 and wave-vector y
  double x0(const double y);  
  // Integrand for for any frequency and wave-vector y
  double integrand(const double y);
  // Integrator object
  const shared_ptr<Integrator1D> itg;
  
public:

  // Constructor for finite temperature calculations
  Idr(const int l_,
      const double x_,
      const double Theta_,
      const double mu_,
      const double yMin_,
      const double yMax_,
      const shared_ptr<Integrator1D> &itg_)
    : l(l_), x(x_), Theta(Theta_),
      mu(mu_), yMin(yMin_), yMax(yMax_),
      itg(itg_) {;};
  // Constructor for zero temperature calculations
  Idr(const double Omega_,
      const double x_)
    : Omega(Omega_), x(x_) {;};
  // Get at finite temperature
  double get();
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
  // Integrator object
  const shared_ptr<Integrator1D> itg;
  // Get integrand
  double integrand(double y);
  
public:

  // Constructor for finite temperature calculations
  SsfHF(const double x_,
	const double Theta_,
	const double mu_,
	const double yMin_,
	const double yMax_,
	const shared_ptr<Integrator1D> &itg_)
    : x(x_), Theta(Theta_), mu(mu_),
      yMin(yMin_), yMax(yMax_), itg(itg_) {;};
  // Constructor for zero temperature calculations
  SsfHF(const double x_) : x(x_) {;};
  // Get at finite temperature
  double get();
  // Get at zero temperature
  double get0();
  
};

// Class for the static local field correction
class Slfc {

private:
  
  // Wave-vector
  const double x = 0;
  // Integration limits
  const double yMin = 0;
  const double yMax = 0;
  // Integrator object
  const shared_ptr<Integrator1D> itg;
  // Integrand
  double integrand(const double y);
  // Static structure factor interpolator
  const shared_ptr<Interpolator> ssfi;
  // Compute static structure factor
  double ssf(double x_);
  
public:

  // Constructor
  Slfc(const double x_,
       const double yMin_,
       const double yMax_,
       const shared_ptr<Integrator1D> &itg_,
       const shared_ptr<Interpolator> &ssfi_)
    : x(x_), yMin(yMin_), yMax(yMax_),
      itg(itg_), ssfi(ssfi_) {;};
  // Get result of integration 
  double get();
  
};

#endif
