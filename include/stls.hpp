#ifndef STLS_HPP
#define STLS_HPP

#include <vector>
#include "input.hpp"
#include "util.hpp"
#include "numerics.hpp"

// -----------------------------------------------------------------
// Solver for the STLS-based schemes
// -----------------------------------------------------------------

class Stls {

protected: 
  
  // Input data
  const StlsInput in;
  // Wave vector grid
  vector<double> wvg;
  // Ideal density response
  vecUtil::Vector2D idr;
  // Static local field correction
  vector<double> slfcOld;
  vector<double> slfc;
  // Static structure factor
  vector<double> ssf;
  // Hartree-Fock static structure factor
  vector<double> ssfHF;
  // Integrator object
  Integrator1D itg;
  // Output verbosity
  const bool verbose;
  // Name of the recovery files
  string recoveryFileName;
  // Flag to write the recovery files
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
  // Write recovery files
  void writeRecovery();
  void readRecovery(vector<double> &wvgFile,
		    vector<double> &slfcFile) const;
  // Check if iet schemes should be used
  bool checkIet() {
   return in.getTheory() == "STLS-HNC" ||
     in.getTheory() == "STLS-IOI" ||
     in.getTheory() == "STLS-LCT";
  };
  
public:

  // Constructors
  Stls(const StlsInput &in_,
       const bool &verbose_,
       const bool &writeFiles_)
    : in(in_), verbose(verbose_), recoveryFileName(EMPTY_STRING),
      writeFiles(writeFiles_), computedChemicalPotential(false),
      useIet(checkIet()) { ; };
  Stls(const StlsInput &in_) : Stls(in_, true, true) { ; };
  // Compute stls scheme
  int compute();
  // Getters
  vector<double> getBf() const { return bf; }
  vecUtil::Vector2D getIdr() const { return idr; }
  string getRecoveryFileName() const { return recoveryFileName; }
  vector<double> getSlfc() const { return slfc; }
  vector<double> getSsf() const { return ssf; }
  vector<double> getSsfHF() const { return ssfHF; }
  vector<double> getWvg() const { return wvg; }
  vector<double> getRdf(const vector<double> &r) const;
  vector<double> getSdr() const;
  double getUInt() const;
  
};

// -----------------------------------------------------------------
// Classes for the ideal density response
// -----------------------------------------------------------------

class Idr {

private:
  
  // Number of matsubara frequency
  const int nl;
  // Wave-vector
  const double x;
  // Degeneracy parameter
  const double Theta;
  // Chemical potential
  const double mu;
  // Integration limits for finite temperature calculations
  const double yMin;
  const double yMax;
  // Idr integrand for frequency = l and wave-vector x
  double integrand(const double y, const int l) const;
  // Idr integrand for frequency = 0 and wave-vector x
  double integrand(const double y) const;  
  // Integrator object
  Integrator1D &itg;
  
public:

  // Constructor
  Idr(const int nl_,
      const double x_,
      const double Theta_,
      const double mu_,
      const double yMin_,
      const double yMax_,
      Integrator1D &itg_)
    : nl(nl_), x(x_), Theta(Theta_),
      mu(mu_), yMin(yMin_), yMax(yMax_),
      itg(itg_) {;};
  // Get result of integration
  vector<double> get() const;
  
};

class IdrGround {

private:
  
  // Frequency
  const double Omega;
  // Wave-vector
  const double x;
  
public:

  // Constructor
  IdrGround(const double Omega_,
	    const double x_) : Omega(Omega_), x(x_) {;};
  // Get real part
  double re0() const;
  // Get imaginary part 
  double im0() const;
  // Get frequency derivative of the real part
  double re0Der() const;
  
};

// -----------------------------------------------------------------
// Classes for the Hartree-Fock static structure factor
// -----------------------------------------------------------------

class SsfHF {

private:
  
  // Wave-vector
  const double x;
  // Degeneracy parameter
  const double Theta;
  // Chemical potential
  const double mu;
  // Integration limits for finite temperature calculations
  const double yMin;
  const double yMax;
  // Integrator object
  Integrator1D &itg;
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
	Integrator1D &itg_) : x(x_), Theta(Theta_), mu(mu_),
			      yMin(yMin_), yMax(yMax_),
			      itg(itg_) {;};
  // Get at any temperature
  double get() const;
  
};

class SsfHFGround {

private:
  
  // Wave-vector
  const double x;
  
public:

  // Constructor for zero temperature calculations
  SsfHFGround(const double x_) : x(x_) {;};
  // Get result
  double get() const;
  
};

// -----------------------------------------------------------------
// Classes for the static structure factor
// -----------------------------------------------------------------

class SsfBase {

protected:
  
  // Wave-vector
  const double x;
  // Degeneracy parameter
  const double Theta;
  // Coupling parameter
  const double rs;
  // Hartree-Fock contribution
  const double ssfHF;
  // Static local field correction
  const double slfc;
  // Constant for unit conversion
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  // Constructor
  SsfBase(const double x_,
	  const double Theta_,
	  const double rs_,
	  const double ssfHF_,
	  const double slfc_) : x(x_), Theta(Theta_), rs(rs_),
				ssfHF(ssfHF_), slfc(slfc_) {;};
  
};

class Ssf : public SsfBase {

protected:

  // Number of Matsubara frequencies
  const int nl;
  // Ideal density response
  const double *idr;
  
public:

  // Constructor
  Ssf(const double x_,
      const double Theta_,
      const double rs_,
      const double ssfHF_,
      const double slfc_,
      const int nl_,
      const double *idr_) : SsfBase(x_, Theta_, rs_, ssfHF_, slfc_),
			    nl(nl_), idr(idr_) {;};
  // Get static structore factor
  double get() const;
 
  
};


class SsfGround : public SsfBase {

private:
  
  // Integration limits for zero temperature calculations
  const double yMin;
  const double yMax;
  // Integrator object
  Integrator1D &itg;
  // Integrand for zero temperature calculations
  double integrand(const double Omega) const ;
  // Plasmon contribution
  double plasmon() const;
  // Dielectric response function
  double drf(const double Omega) const;
  // Frequency derivative of the dielectric response function
  double drfDer(const double Omega) const;
  
public:

  // Constructor for zero temperature calculations
  SsfGround(const double x_,
	    const double rs_,
	    const double ssfHF_,
	    const double slfc_,
	    const double yMin_,
	    const double yMax_,
	    Integrator1D &itg_) : SsfBase(x_, 0, rs_, ssfHF_, slfc_),
				  yMin(yMin_), yMax(yMax_), itg(itg_) {;};
  // Get result of integration
  double get() const;
 
  
};

// -----------------------------------------------------------------
// Classes for the static local field correction
// -----------------------------------------------------------------

class SlfcBase {

protected:
  
  // Wave-vector
  const double x;
  // Integration limits
  const double yMin;
  const double yMax;
  // Static structure factor interpolator
  const Interpolator1D &ssfi;
  // Compute static structure factor
  double ssf(double x_) const;
  // Constructor
  SlfcBase(const double x_,
	   const double yMin_,
	   const double yMax_,
	   const Interpolator1D &ssfi_) : x(x_), yMin(yMin_),
					  yMax(yMax_), ssfi(ssfi_) {;};
  
};

class Slfc : public SlfcBase {

private:
  
  // Integrator object
  Integrator1D &itg;
  // Integrand
  double integrand(const double y) const;
  
public:

  // Constructor
  Slfc(const double x_,
       const double yMin_,
       const double yMax_,
       const Interpolator1D &ssfi_,
       Integrator1D &itg_) : SlfcBase(x_, yMin_, yMax_, ssfi_),
			     itg(itg_) { ; };
  // Get result of integration 
  double get() const;
  
};


class SlfcIet : public SlfcBase {

private:

  // Integrator object
  Integrator2D &itg;
  // Grid for 2D integration
  const vector<double> &itgGrid;
  // Integrands
  double integrand1(const double y) const;
  double integrand2(const double w) const;
  // Static local field correction interpolator
  const Interpolator1D &slfci;
  // Bridge function interpolator
  const Interpolator1D &bfi;
  // Compute static local field correction
  double slfc(double x_) const;
  // Compute bridge function
  double bf(double x_) const;
  
public:

  // Constructor
  SlfcIet(const double x_,
	  const double yMin_,
	  const double yMax_,
	  const Interpolator1D &ssfi_,
	  const Interpolator1D &slfci_,
	  const Interpolator1D &bfi_,
	  const vector<double> &itgGrid_,
	  Integrator2D &itg_)
    : SlfcBase(x_, yMin_, yMax_, ssfi_), itg(itg_),
      itgGrid(itgGrid_), slfci(slfci_), bfi(bfi_) {;};
  // Get result of integration 
  double get() const;

};


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
  Integrator1DFourier &itg;
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
  BridgeFunction(const string theory_,
		 const string mapping_,
		 const double rs_,
		 const double Theta_,
		 const double x_,
		 Integrator1DFourier &itg_) : theory(theory_),
					      mapping(mapping_),
					      rs(rs_), Theta(Theta_),
					      x(x_), itg(itg_) {;};
  // Get result of the integration
  double get() const;
  
};

#endif
