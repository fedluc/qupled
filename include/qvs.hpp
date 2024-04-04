#ifndef QVS_HPP
#define QVS_HPP

#include <limits>
#include <map>
#include <functional>
#include "numerics.hpp"
#include "vsbase.hpp"
#include "qstls.hpp"

// Forward declarations
class QVSStlsInput;

// -----------------------------------------------------------------
// Class to handle simultaneous state point calculations
// -----------------------------------------------------------------

class QStlsCSR : public CSR<vecUtil::Vector2D, Qstls, QVSStlsInput> {

  friend class QStructProp;

private:
  
  // Compute auxiliary density response
  void computeAdrStls();
  void computeAdr();
  // Helper methods to compute the derivatives
  double getDerivative(const vecUtil::Vector2D& f,
		       const int &l,
		       const size_t& idx,
		       const Derivative& type);

public:
  
  // Constructor
  QStlsCSR(const QVSStlsInput& in_) : CSR(in_, Qstls(in_)) { ; }
  // Compute Q
  double getQAdder() const;
};

// -----------------------------------------------------------------
// Class to handle the Q-adder in the free parameter expression
// -----------------------------------------------------------------

class QAdder {

private:
  
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  // --Constructor variables--
  // Coupling parameter
  const double rs;
  // Degeneracy parameter
  const double Theta;
  // Chemical potential
  const double mu;
  // Integration limits
  const std::pair<double, double> limits;
  // Grid for 2D integration
  const std::vector<double> &itgGrid;
  
  // Integrator1D 
  Integrator1D &itg1;
  // Integrator2D 
   Integrator2DSingular &itg2;
  // Interpolator 1D class instance
  const Interpolator1D &interp;
  
  // SSF interpolation
  double ssf(const double& y) const;
  // Integrands 
  double integrandDenominator(const double q) const;
  double integrandNumerator1(const double q) const;
  double integrandNumerator2(const double w) const;
  // Get Integral denominator
  void getIntDenominator(double &res) const;
  
  
public:

  // Constructor
  QAdder(const double& rs_,
	 const double& Theta_,
	 const double& mu_,
	 double limitMin, double limitMax,
	 const std::vector<double>& itgGrid_,
	 Integrator1D &itg1_,
	 Integrator2DSingular &itg2_,
	 const Interpolator1D &interp_) :  rs(rs_), Theta(Theta_), mu(mu_), 
             limits(limitMin, limitMax), itgGrid(itgGrid_), itg1(itg1_), 
             itg2(itg2_), interp(interp_) { ; }
  
  // Get Q-adder
  double get() const;
  
};

// -----------------------------------------------------------------
// Class to handle the state point derivatives
// -----------------------------------------------------------------

class QStructProp : public StructPropBase<QStlsCSR, QVSStlsInput> {

private:

  using StructPropBase = StructPropBase<QStlsCSR, QVSStlsInput>;
  // Perform iterations to compute structural properties
  void doIterations();
  
  
public:

  // Compute structural properties
  int compute();
  // Constructor
  QStructProp(const QVSStlsInput &in_) : StructPropBase(in_) { ; }
  // Get Q term
  std::vector<double> getQ() const;  
  
};


// -----------------------------------------------------------------
// Class to handle the thermodynamic properties and derivatives
// -----------------------------------------------------------------

class QThermoProp : public ThermoPropBase<QStructProp, QVSStlsInput> {

private:

  using ThermoPropBase = ThermoPropBase<QStructProp, QVSStlsInput>;

public:

  // Constructors
  QThermoProp(const QVSStlsInput& in) : ThermoPropBase(in) { ; } 
  QThermoProp(const QVSStlsInput& in,
	      const QThermoProp& other) : ThermoPropBase(in, other) { ; }
  QThermoProp(const ThermoPropBase& other) : ThermoPropBase(other) { ; }
  // Get internal energy and internal energy derivatives
  std::vector<double> getQData() const;
  
};

// -----------------------------------------------------------------
// Solver for the qVS-STLS scheme
// -----------------------------------------------------------------

class QVSStls : public VSBase<QThermoProp, Rpa, QVSStlsInput> {

private:

  using VSBase = VSBase<QThermoProp, Rpa, QVSStlsInput>;
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the qvsstls-scheme
  void updateSolution();
  // Auxiliary density response
  vecUtil::Vector2D adr;
  
public:
  
  // Constructor from initial data
  QVSStls(const QVSStlsInput &in_) : VSBase(in_) { ; }
  // Constructor for recursive calculations
  QVSStls(const QVSStlsInput &in_,
	  const QThermoProp& thermoProp_) : VSBase(in_, thermoProp_) { ; }

  // Getters
  vecUtil::Vector2D getAdr() const { return adr; }
  
};

#endif
