#ifndef STLS_HPP
#define STLS_HPP

#include <vector>
#include "input.hpp"
#include "numerics.hpp"

using namespace std;
using namespace vecUtil;

class Stls {

protected: 

  // Wave vector grid
  vector<double> wvg;
  // Ideal density response
  Vector2D<double> idr;
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
  // Flag to write the outputfiles
  const bool writeFiles;
  // Chemical potential
  double mu;
  bool computedChemicalPotential;
  // iet schemes
  bool useIet;
  // Bridge function (for iet schemes)
  vector<double> bf;
  // Constant for unit conversion
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  // Initialization
  void init();
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
  void computeSlfcStls();
  void computeSlfcIet();
  // Compute bridge function
  void computeBf();
  // Iterations to solve the stls scheme
  void doIterations();
  void initialGuess();
  double computeError();
  void updateSolution();
  // Write output files
  void writeOutput() const;
  void writeSsf() const;
  void writeSsfHF() const;
  void writeSlfc() const;
  void writeSdr() const;
  void writeIdr() const;
  void writeUInt() const;
  void writeRdf() const;
  void writeBf() const;
  // Restart files
  void writeRestart() const;
  void readRestart(vector<double> &wvgFile,
		   vector<double> &slfcFile) const;
  // Check if iet schemes should be used
  void checkIet() { useIet = in.getTheory() == "STLS-HNC" ||
      in.getTheory() == "STLS-IOI" ||
      in.getTheory() == "STLS-LCT";}
  
public:

  // Constructors
  Stls(const Input in_)
    : in(in_), verbose(true),
      writeFiles(true), computedChemicalPotential(false) {checkIet();};
  Stls(const Input in_, const bool verbose_, const bool writeFiles_)
    : in(in_), verbose(verbose_),
      writeFiles(writeFiles_), computedChemicalPotential(false) {checkIet();};
  // Compute stls scheme
  void compute();
  // Getters
  void getSsf(vector<double> &ssf_);

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
  double integrand(const double y, const int l) const;
  // Idr integrand for frequency = 0 and wave-vector x
  double integrand(const double y) const;  
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
  vector<double> get() const;
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
  double integrand(double y) const;
  // Get at zero temperature
  double get0() const;
  
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
  double get() const;
  
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
  // Auxiliary density response
  const vector<double> adr;
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
  double integrand(const double Omega) const ;
  // Plasmon contribution
  double plasmon() const;
  // Dielectric response function
  double drf(const double Omega) const;
  // Frequency derivative of the dielectric response function
  double drfDer(const double Omega) const;
  // Get at zero temperature
  double get0() const;
  // Get for quantum schemes
  double getQuantum() const;
  
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
  // Constructor for finite temperature calculations for quantum schemes
  Ssf(const double x_,
      const double Theta_,
      const double rs_,
      const double ssfHF_,
      const vector<double> &idr_,
      const vector<double> &adr_)
    : x(x_), Theta(Theta_), rs(rs_),
      ssfHF(ssfHF_), idr(idr_), adr(adr_) {;};
  // Get at any temperature
  double get() const;
 
  
};

// Class for the static local field correction
class Slfc {

protected:
  
  // Wave-vector
  const double x = 0;
  // Integration limits
  const double yMin = 0;
  const double yMax = 0;

private:
  
  // Integrator object
  const shared_ptr<Integrator1D> itg;
  // Integrand
  double integrand(const double y) const;

protected:
  
  // Static structure factor interpolator
  const shared_ptr<Interpolator> ssfi;
  // Compute static structure factor
  double ssf(double x_) const;
  
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
  double get() const;
  
};


// Class for the static local field correction for the iet schemes
class SlfcIet : public Slfc {

private:

  // Integrator object
  const shared_ptr<Integrator2D> itg2;
  // Integrands
  double integrand1(const double y) const;
  double integrand2(const double w) const;
  // Static local field correction interpolator
  const shared_ptr<Interpolator> slfci;
  // Bridge function interpolator
  const shared_ptr<Interpolator> bfi;
  // Compute static local field correction
  double slfc(double x_) const;
  // Compute bridge function
  double bf(double x_) const;
  
public:

  // Constructor
  SlfcIet(const double x_,
	  const double yMin_,
	  const double yMax_,
	  const shared_ptr<Integrator2D> &itg2_,
	  const shared_ptr<Interpolator> &ssfi_,
	  const shared_ptr<Interpolator> &slfci_,
	  const shared_ptr<Interpolator> &bfi_)
    : Slfc(x_, yMin_, yMax_, NULL, ssfi_),
      itg2(itg2_), slfci(slfci_), bfi(bfi_) {;};
  // Get result of integration 
  double get() const;

};


// Class for the bridge function
class BridgeFunction {

private:

  // Theory to be solved
  const string theory;
  // Iet mapping
  const string mapping;
  // Coupling parameter
  const double rs;
  // Degeneracy parameter
  const double Theta;
  // Wave vector
  const double x;
  // Integrator object
  const shared_ptr<Integrator1DFourier> itg;
  // Constant for unit conversion
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  // Hypernetted-chain bridge function
  double hnc() const ;
  // Ichimaru bridge function
  double ioi() const;
  // Lucco Castello and Tolias bridge function
  double lct() const;
  double lctIntegrand(const double r, const double Gamma) const;
  // Coupling parameter to compute the bridge function
  double couplingParameter() const;

public:

  // Constructor
  BridgeFunction(const string theory_, const string mapping_,
		 const double rs_, const double Theta_,
		 const double x_, const shared_ptr<Integrator1DFourier> &itg_)
    : theory(theory_), mapping(mapping_), rs(rs_), Theta(Theta_),
      x(x_), itg(itg_) {;};
  // getter
  double get() const;
  
};

#endif
