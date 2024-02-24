#ifndef QVS_HPP
#define QVS_HPP

#include <limits>
#include <map>
#include <functional>
#include "numerics.hpp"
#include "vsstls.hpp"
#include "qstls.hpp"


// // -----------------------------------------------------------------
// // Class to handle simultaneous state point calculations
// // -----------------------------------------------------------------

class qStlsCSR : public Qstls, public CSR<vecUtil::Vector2D> {

  friend class qStructProp;

private:

  // Input data
  const VSStlsInput in;
  // Compute auxiliary density response
  void computeAdrStls();
  void computeAdr();
  // Compute the internal energy
  double getInternalEnergy() const;
  // Compute the free energy integrand
  double getFreeEnergyIntegrand() const;
  // Compute Q
  double getQ() const;
  // Convert VSStlsInput to QstlsInput
  QstlsInput VStoQStlsInput(const VSStlsInput& in) const;

public:
  
  // Constructor
  qStlsCSR(const VSStlsInput& in_) : Qstls(VStoQStlsInput(in_)),
				     in(in_) { ; }
  
  
};


// -----------------------------------------------------------------
// Class to handle the state point derivatives
// -----------------------------------------------------------------

class qStructProp : public StructProp {

private:

  // Vector containing NPOINTS state points to be solved simultaneously
  std::vector<qStlsCSR> Adr;
  // Flag marking whether the initialization for the stls data is done
  bool adrIsInitialized;
  // Perform iterations to compute structural properties
  void doIterations();
  // Generic getter function to return vector data
  const std::vector<double>& getBase(std::function<double(const qStlsCSR&)> f) const;
  
public:

  // Constructor
  qStructProp(const VSStlsInput &in_);
  // Compute structural properties
  int compute();
  // Set free parameter
  void setAlpha(const double& alpha);
  // Get coupling parameters for all the state points
  std::vector<double> getCouplingParameters() const;
  // Get degeneracy parameters for all the state points
  std::vector<double> getDegeneracyParameters() const;
  // Get internal energy for all the state points
  std::vector<double> getInternalEnergy() const;
  // Get free energy integrand for all the state points
  std::vector<double> getFreeEnergyIntegrand() const;
  // Get Q term
  std::vector<double> getQ() const;
  const qStlsCSR& getAdrStls(const Idx& idx) const { return Adr[idx]; }
  
  
};


// -----------------------------------------------------------------
// Class to handle the thermodynamic properties and derivatives
// -----------------------------------------------------------------

class qThermoProp : public ThermoProp {

private:

  // // Compute the free energy
  // double computeQ(const SIdx iStruct,
  // 		  const bool normalize) const ;
  qStructProp qstructProp;
  
public:

  // Constructors
  qThermoProp(const VSStlsInput& in) : ThermoProp(in), qstructProp(in) { ; } 
  qThermoProp(const VSStlsInput& in,
	      const qThermoProp& other) : ThermoProp(in, other), qstructProp(in) { ; } 
  void compute(const VSStlsInput& in);
  // Get internal energy and internal energy derivatives
  std::vector<double> getQData() const;
  // Set the value of the free parameter in the structural properties
  void setAlpha(const double& alpha);
  const qStlsCSR& getStructProp();
  
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

class qVSStls : public VSStls {

private:

  qThermoProp qthermoProp;
  double computeAlpha();
  double alphaDifference(const double& alphaTmp);
  void updateSolution();
  // Auxiliary density response
  vecUtil::Vector2D adr;
  
public:
  
  
  // Constructor from initial data
  qVSStls(const VSStlsInput &in_) : VSStls(in_),
				    qthermoProp(in_) { ; }
  // Constructor for recursive calculations
  qVSStls(const VSStlsInput &in_,
	  const qThermoProp& qthermoProp_) : VSStls(in_, qthermoProp_),
					     qthermoProp(in_, qthermoProp_) { ; }
  // Getters
  const qThermoProp& getThermoProp() const { return qthermoProp; }
  
};

#endif
