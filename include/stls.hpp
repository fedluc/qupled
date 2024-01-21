#ifndef STLS_HPP
#define STLS_HPP

#include <vector>
#include "input.hpp"
#include "util.hpp"
#include "numerics.hpp"
#include "rpa.hpp"

// -----------------------------------------------------------------
// Solver for the STLS-based schemes
// -----------------------------------------------------------------

#define EMPTY_STRING ""

class StlsBase {

protected:

  // Constant for unit conversion
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
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
  // Bridge function (for iet schemes)
  vector<double> bf;
  // iet schemes
  bool useIet;
  // Chemical potential
  double mu;
  bool computedChemicalPotential;
  // Name of the recovery files
  string recoveryFileName;
  

public:

  // Constructor
  StlsBase(const StlsInput &in_) : in(in_), useIet(false),
				   computedChemicalPotential(false),
				   recoveryFileName(EMPTY_STRING) { ; };
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

class Stls : public StlsBase {

protected: 
  
  // Integrator object
  Integrator1D itg;
  // Output verbosity
  const bool verbose;
  // Flag to write the recovery files
  const bool writeFiles;
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
  void checkIet() {
   useIet = in.getTheory() == "STLS-HNC"
     || in.getTheory() == "STLS-IOI"
     || in.getTheory() == "STLS-LCT";
  };
  
public:

  // Constructors
  Stls(const StlsInput &in_,
       const bool &verbose_,
       const bool &writeFiles_)
    : StlsBase(in_), itg(in_.getIntError()), verbose(verbose_),
      writeFiles(writeFiles_) { checkIet(); };
  Stls(const StlsInput &in_) : Stls(in_, true, true) { ; };
  // Compute stls scheme
  int compute();

  
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
