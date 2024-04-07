#ifndef RPA_HPP
#define RPA_HPP

#include <vector>

// Forward declarations
namespace vecUtil {
  class Vector2D;
}
class RpaInput;
class IntegratorCQUAD;

// -----------------------------------------------------------------
// Solver for the Random-Phase approximation scheme
// -----------------------------------------------------------------

class Rpa {

protected:
  
  // Constant for unit conversion
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  // Input data
  const RpaInput in;
  // Output verbosity
  const bool verbose;
  // Name of the recovery files
  std::string recoveryFileName;
  // Integrator
  IntegratorCQUAD itg;
  // Wave vector grid
  std::vector<double> wvg;
  // Ideal density response
  vecUtil::Vector2D idr;
  // Static local field correction
  std::vector<double> slfc;
  // Static structure factor
  std::vector<double> ssf;
  // Hartree-Fock static structure factor
  std::vector<double> ssfHF;
  // Chemical potential
  double mu;
  // Initialize basic properties
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

public:

  // Constructor
  Rpa(const RpaInput &in_,
      const bool verbose_);
  Rpa(const RpaInput& in_) : Rpa(in_, true) { ; }
  // Compute the scheme
  int compute();
  // Getters
  vecUtil::Vector2D getIdr() const { return idr; }
  std::vector<double> getSlfc() const { return slfc; }
  std::vector<double> getSsf() const { return ssf; }
  std::vector<double> getSsfHF() const { return ssfHF; }
  std::vector<double> getWvg() const { return wvg; }
  std::vector<double> getRdf(const std::vector<double>& r) const;
  std::vector<double> getSdr() const;
  double getUInt() const;
  std::string getRecoveryFileName() const { return recoveryFileName; }
  
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
  double integrand(const double& y,
		   const int& l) const;
  // Idr integrand for frequency = 0 and wave-vector x
  double integrand(const double& y) const;  
  // Integrator object
  IntegratorCQUAD &itg;
  
public:

  // Constructor
  Idr(const int nl_,
      const double& x_,
      const double& Theta_,
      const double& mu_,
      const double& yMin_,
      const double& yMax_,
      IntegratorCQUAD &itg_) : nl(nl_), x(x_), Theta(Theta_),
			       mu(mu_), yMin(yMin_), yMax(yMax_),
			       itg(itg_) {;};
  // Get result of integration
  std::vector<double> get() const;
  
};

class IdrGround {

private:
  
  // Frequency
  const double Omega;
  // Wave-vector
  const double x;
  
public:

  // Constructor
  IdrGround(const double& Omega_,
	    const double& x_) : Omega(Omega_), x(x_) {;};
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
  IntegratorCQUAD &itg;
  // Get integrand
  double integrand(const double& y) const;
  // Get at zero temperature
  double get0() const;
  
public:

  // Constructor for finite temperature calculations
  SsfHF(const double& x_,
	const double& Theta_,
	const double& mu_,
	const double& yMin_,
	const double& yMax_,
	IntegratorCQUAD &itg_) : x(x_), Theta(Theta_), mu(mu_),
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
  SsfHFGround(const double& x_) : x(x_) {;};
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
  SsfBase(const double& x_,
	  const double& Theta_,
	  const double& rs_,
	  const double& ssfHF_,
	  const double& slfc_) : x(x_), Theta(Theta_), rs(rs_),
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
  Ssf(const double& x_,
      const double& Theta_,
      const double& rs_,
      const double& ssfHF_,
      const double& slfc_,
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
  IntegratorCQUAD &itg;
  // Integrand for zero temperature calculations
  double integrand(const double& Omega) const ;
  // Plasmon contribution
  double plasmon() const;
  // Dielectric response function
  double drf(const double& Omega) const;
  // Frequency derivative of the dielectric response function
  double drfDer(const double& Omega) const;
  
public:

  // Constructor for zero temperature calculations
  SsfGround(const double& x_,
	    const double& rs_,
	    const double& ssfHF_,
	    const double& slfc_,
	    const double& yMin_,
	    const double& yMax_,
	    IntegratorCQUAD &itg_) : SsfBase(x_, 0, rs_, ssfHF_, slfc_),
				     yMin(yMin_), yMax(yMax_), itg(itg_) {;};
  // Get result of integration
  double get() const;
 
  
};

#endif
