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
  // Compute auxiliary density response
  void computeAdrStls();
  void computeAdr();
  // Compute Q
  double getQ() const;

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

// class Q : public AdrFixed {
//     // Get integration result
//     void get(std::vector<double> &wvg,
// 	   vecUtil::Vector3D &res) const;
//     double integranddenom(const double q,
// 		    const double l) const;
//     double integrandnum(const double t,
// 		    const double y,
// 		    const double l) const;

// };

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
