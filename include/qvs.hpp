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

class qStlsCSR : public CSR<vecUtil::Vector2D, Qstls, QVSStlsInput> {

  friend class qStructProp;

private:
  
  // Compute auxiliary density response
  void computeAdrStls();
  void computeAdr();
  // Compute Q
  double getQ() const;
  // Helper methods to compute the derivatives
  double getDerivative(const vecUtil::Vector2D& f,
		       const int &l,
		       const size_t& idx,
		       const Derivative& type);

public:
  
  // Constructor
  qStlsCSR(const QVSStlsInput& in_) : CSR(in_, Qstls(in_)) { ; }
  
};

// // -----------------------------------------------------------------
// // Class to handle the Q-adder in the free parameter expression
// // -----------------------------------------------------------------

// // COMMENT: Give a different name to this class, single letter names for
// // classes can be problematic when maintining the code (e.g. if you try
// // to search for Q to see where Q is used you will get a lot of
// // irrelevant results). Single letter names for variables are ok if the
// // scope of the variable is short

// class Q {

// private:

//   // COMMENT: You don't need this object if you already pass Theta, rs and mu
  
//   const VSStlsInput in;
  
//   // COMMENT: You call these constructor variables but you haven't defined
//   // a constructor that sets them. Moreover, constant members can only be set in
//   // the initialization list of the constructor (see my comment in the public
//   // section)
  
//   // --Constructor variables--
//   // Degeneracy parameter
//   const double Theta;
//   // Coupling parameter
//   const double rs;
//   // Chemical potential
//   const double mu;
  
//   // COMMENT: Instead of computing a fixed part use an Integrator2D
//   // object to compute the numerator. If you want an example on how to use
//   // it take a look at the SlfcIet class. In that class we are solving
//   // a double integral to compute the slfc, the double integral that is being
//   // solved can be found here:
//   // https://pubs.aip.org/aip/jcp/article/155/13/134115/353165/Integral-equation-theory-based-dielectric-scheme.
  
//   // Vector value for Numerator fixed part
//   std::vector<double> NumPart;
//   // Value for Numerator
//   double Numerator;
//   // Value for Denominator
//   double Denominator;
  
//   // COMMENT: Here Integrator1D and Interpolator1D are reference members,
//   // make sure you initialize them in a constructor or the code won't compile.
//   // If you want to make your life easier make them regular members by
//   // removing the & symbol. In other parts of the code I use references
//   // for these variables because I want to avoid repeated memory allocations
//   // when the integrals are compute in a for loop,
//   // but here you compute the integral only once per iteration so
//   // memory allocations shouldn't be an issue.
  
//   // Integrator object instance
//   Integrator1D &itg;
//   // Interpolator 1D class instance
//   Interpolator1D &interp;
  
//   // COMMENT: Make this a constexpr
  
//   const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);

//   // Integrands 
//   double integrandDenominator(const double q) const;
//   double integrandNumerator(const double t,
// 		    const double y) const;
//   // Get Integrals
//   void getIntNumerator(const vector<double> wvg,
//         std::vector<double> &res) const;
//   void getIntDenominator(const vector<double> wvg,
//         double &res) const;
//   void getTotalNumerator(const vector<double> wvg,
//       std::vector<double> &res,
//       std::vector<double> ssf,
//       double &Num) const;

  
// public:

//   // COMMENT: Add a constructor where you initialize Theta, rs, mu.
//   // Something along this line:
//   // Q(const double& rs_,
//   //   const double& Theta_,
//   //   const double& mu_) : rs(rs_), Theta(Theta_), mu(mu_) { ; }

//   // COMMENT: You can rename this compute or get (check what I use
//   // for Slfc or similar). Also, since rs and Theta are defined in the
//   // constructor then you don't need to pass them here 
//   double computeQ(const vector<double> &wvg, 
//                   const std::vector<double> ssf, 
//                   const double rs, 
//                   const double Theta);
  
// };

// -----------------------------------------------------------------
// Class to handle the state point derivatives
// -----------------------------------------------------------------

class qStructProp : public StructPropBase<qStlsCSR, QVSStlsInput> {

private:

  using StructPropBase = StructPropBase<qStlsCSR, QVSStlsInput>;
  // Perform iterations to compute structural properties
  void doIterations();
  
public:

  // Constructor
  qStructProp(const QVSStlsInput &in_) : StructPropBase(in_) { ; }
  // Get Q term
  std::vector<double> getQ() const;  
  
};


// -----------------------------------------------------------------
// Class to handle the thermodynamic properties and derivatives
// -----------------------------------------------------------------

class qThermoProp : public ThermoPropBase<qStructProp, QVSStlsInput> {

private:

  using ThermoPropBase = ThermoPropBase<qStructProp, QVSStlsInput>;

public:

  // Constructors
  qThermoProp(const QVSStlsInput& in) : ThermoPropBase(in) { ; } 
  qThermoProp(const QVSStlsInput& in,
	      const qThermoProp& other) : ThermoPropBase(in, other) { ; }
  qThermoProp(const ThermoPropBase& other) : ThermoPropBase(other) { ; }
  // Get internal energy and internal energy derivatives
  std::vector<double> getQData() const;
  
};

// -----------------------------------------------------------------
// Solver for the qVS-STLS scheme
// -----------------------------------------------------------------

class qVSStls : public VSBase<qThermoProp, Rpa, QVSStlsInput> {

private:

  using VSBase = VSBase<qThermoProp, Rpa, QVSStlsInput>;
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the qvsstls-scheme
  void updateSolution();
  // Auxiliary density response
  vecUtil::Vector2D adr;
  
public:
  
  // Constructor from initial data
  qVSStls(const QVSStlsInput &in_) : VSBase(in_) { ; }
  // Constructor for recursive calculations
  qVSStls(const QVSStlsInput &in_,
	  const qThermoProp& thermoProp_) : VSBase(in_, thermoProp_) { ; }
  
  
};

#endif
