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

class Stls : public Rpa {

protected: 

  // Input parameters
  StlsInput in;
  // Flag to write the recovery files
  const bool writeFiles;
  // iet schemes
  bool useIet;
  // Name of the recovery files
  string recoveryFileName;
  // Static local field correction to use during the iterations
  vector<double> slfcOld;
  // Bridge function (for iet schemes)
  vector<double> bf;
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
  
public:

  // Constructors
  Stls(const StlsInput &in_,
       const bool verbose_,
       const bool writeFiles_);
  Stls(const StlsInput &in_) : Stls(in_, true, true) { ; };
  // Compute stls scheme
  int compute();
  // Getters
  vector<double> getBf() const { return bf; }
  string getRecoveryFileName() const { return recoveryFileName; }
  
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
