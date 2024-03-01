#ifndef QVS_HPP
#define QVS_HPP

#include <limits>
#include <map>
#include <functional>
#include "numerics.hpp"
#include "vsbase.hpp"
#include "qstls.hpp"


// -----------------------------------------------------------------
// Class to handle simultaneous state point calculations
// -----------------------------------------------------------------

class qStlsCSR : public CSR<vecUtil::Vector2D, Qstls> {

  friend class qStructProp;

private:
  
  // Input data
  const VSStlsInput in;
  Q QInstance;
  // Compute auxiliary density response
  void computeAdrStls();
  void computeAdr();
  // Compute Q
  double getQ() const;
  // Convert VSStlsInput to QstlsInput
  QstlsInput VStoQStlsInput(const VSStlsInput& in) const;
  double getqDerivative(const vecUtil::Vector2D& f,
            const int &l,
            const size_t& idx,
            const Derivative& type);

public:
  
  // Constructor
  qStlsCSR(const VSStlsInput& in_) : CSR(in_, Qstls(QstlsInput(in_))) { ; }
  
};


// -----------------------------------------------------------------
// Class to handle the state point derivatives
// -----------------------------------------------------------------

class qStructProp : public StructPropBase<qStlsCSR> {

private:

  using StructPropBase = StructPropBase<qStlsCSR>;
  // Perform iterations to compute structural properties
  void doIterations();
  
public:

  // Constructor
  qStructProp(const VSStlsInput &in_) : StructPropBase(in_) { ; }
  // Get Q term
  std::vector<double> getQ() const;  
  
};


// -----------------------------------------------------------------
// Class to handle the thermodynamic properties and derivatives
// -----------------------------------------------------------------

class qThermoProp : public ThermoPropBase<qStructProp> {

private:

  using ThermoPropBase = ThermoPropBase<qStructProp>;
  // // Compute the free energy
  // double computeQ(const SIdx iStruct,
  // 		  const bool normalize) const ;

public:

  // Constructors
  qThermoProp(const VSStlsInput& in) : ThermoPropBase(in) { ; } 
  qThermoProp(const VSStlsInput& in,
	      const qThermoProp& other) : ThermoPropBase(in, other) { ; }
  qThermoProp(const ThermoPropBase& other) : ThermoPropBase(other) { ; }
  // Get internal energy and internal energy derivatives
  std::vector<double> getQData() const;
  
};

class Q {

private:

  const VSStlsInput in;
  // --Constructor variables--
  // Degeneracy parameter
  const double Theta;
  // Coupling parameter
  const double rs;
  // Chemical potential
  const double mu;
  // Vector value for Numerator fixed part
  std::vector<double> NumPart;
  // Value for Numerator
  double Numerator;
  // Value for Denominator
  double Denominator;
  // Integrator object instance
  Integrator1D &itg;
  // Interpolator 1D class instance
  Interpolator1D &interp;
  const double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);

  // Integrands 
  double integrandDenominator(const double q) const;
  double integrandNumerator(const double t,
		    const double y) const;
  // Get Integrals
  void getIntNumerator(const vector<double> wvg,
        std::vector<double> &res) const;
  void getIntDenominator(const vector<double> wvg,
        double &res) const;
  void getTotalNumerator(const vector<double> wvg,
      std::vector<double> &res,
      std::vector<double> ssf,
      double &Num) const;

  
public:
  
  double computeQ(const vector<double> &wvg, 
                  const std::vector<double> ssf, 
                  const double rs, 
                  const double Theta);
  
};


// -----------------------------------------------------------------
// Solver for the qVS-STLS scheme
// -----------------------------------------------------------------

class qVSStls : public VSBase<qThermoProp, Rpa> {

private:

  using VSBase = VSBase<qThermoProp, Rpa>;
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the qvsstls-scheme
  void updateSolution();
  // Auxiliary density response
  vecUtil::Vector2D adr;
  
public:
  
  // Constructor from initial data
  qVSStls(const VSStlsInput &in_) : VSBase(in_) { ; }
  // Constructor for recursive calculations
  qVSStls(const VSStlsInput &in_,
	  const qThermoProp& thermoProp_) : VSBase(in_, thermoProp_) { ; }
  
  
};

#endif
