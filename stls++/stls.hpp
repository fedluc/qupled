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
  // Integrator object
  const shared_ptr<Integrator1D> itg = make_shared<Integrator1D>();
  // Input data
  const Input in;
  // Output verbosity
  const bool verbose;
  // Chemical potential
  double mu;
  bool computedChemicalPotential;
  // Constant for unit conversion
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  // Construct wave vector grid
  void buildWvGrid();
  // Compute chemical potential
  void computeChemicalPotential();
  // Compute the ideal density response
  void computeIdr();
  // Compute Hartree-Fock static structure factor
  void computeSsfHF();
  void computeSsfHFFinite();
  void computeSsfHFGround();
  // Compute static structure factor at finite temperature
  void computeSsf();
  void computeSsfFinite();
  void computeSsfGround();
  // Compute static local field correction
  void computeSlfc();
  // Iterations to solve the stls scheme
  void doIterations();
  void initialGuess();
  double computeError();
  void updateSolution();
  // Write output files
  void writeOutput();
  void writeSsf();
  void writeSsfHF();
  void writeSlfc();
  void writeSdr();
  void writeIdr();
  void writeRdf();
  void writeUInt();

public:

  // Constructors
  Stls(const Input in_)
    : in(in_), verbose(true), computedChemicalPotential(false) {;};
  Stls(const Input in_, const bool verbose_)
    : in(in_), verbose(verbose_), computedChemicalPotential(false) {;};
  // Compute stls scheme
  void compute(); 

};

// Class for the ideal density response calculation
class Idr {

private:
  
  // Number of matsubara frequency
  const int nl = 0;
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
  double integrand(const double y, const int l);
  // Idr integrand for frequency = 0 and wave-vector x
  double integrand(const double y);  
  // Integrator object
  const shared_ptr<Integrator1D> itg;
  
public:

  // Constructor for finite temperature calculations
  Idr(const int nl_,
      const double x_,
      const double Theta_,
      const double mu_,
      const double yMin_,
      const double yMax_,
      const shared_ptr<Integrator1D> &itg_)
    : nl(nl_), x(x_), Theta(Theta_),
      mu(mu_), yMin(yMin_), yMax(yMax_),
      itg(itg_) {;};
  // Constructor for zero temperature calculations
  Idr(const double Omega_,
      const double x_)
    : Omega(Omega_), x(x_) {;};
  // Get at finite temperature
  vector<double> get();
  // Get real part at zero temperature
  double re0() const;
  // Get imaginary part at zero temperature
  double im0() const;
  // Get frequency derivative of the real part at zero temperature
  double re0Der() const;
  
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
  // Get at zero temperature
  double get0();
  
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
  // Get at any temperature
  double get();
  
};


// Class for the static structure factor
class Ssf {

private:
  
  // Wave-vector
  const double x = 0;
  // Degeneracy parameter
  const double Theta = 0;
  // Coupling parameter
  const double rs = 0;
  // Hartree-Fock contribution
  const double ssfHF = 0;
  // Ideal density response
  const vector<double> idr;
  // Static local field correction
  const double slfc = 0;
  // Constant for unit conversion
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  // Integration limits for zero temperature calculations
  const double yMin = 0;
  const double yMax = 0;
  // Integrator object
  const shared_ptr<Integrator1D> itg;
  // Integrand for zero temperature calculations
  double integrand(const double Omega);
  // Plasmon contribution
  double plasmon();
  // Dielectric response function
  double drf(const double Omega);
  // Frequency derivative of the dielectric response function
  double drfDer(const double Omega);
  // Get at zero temperature
  double get0();
  
public:

  // Constructor for finite temperature calculations
  Ssf(const double x_,
      const double Theta_,
      const double rs_,
      const double ssfHF_,
      const vector<double> &idr_,
      const double slfc_)
    : x(x_), Theta(Theta_), rs(rs_),
      ssfHF(ssfHF_), idr(idr_), slfc(slfc_) {;};
  // Constructor for zero temperature calculations
  Ssf(const double x_,
      const double rs_,
      const double ssfHF_,
      const double slfc_,
      const double yMin_,
      const double yMax_,
      const shared_ptr<Integrator1D> &itg_)
    : x(x_), rs(rs_), ssfHF(ssfHF_),
      slfc(slfc_), yMin(yMin_), yMax(yMax_),
      itg(itg_) {;};
  // Get at any temperature
  double get();
 
  
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
