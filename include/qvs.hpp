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
  double getqDerivative(const Vector2D& f,
            const int &l,
            const size_t& idx,
            const Derivative& type);

public:
  
  // Constructor
  qStlsCSR(const VSStlsInput& in_) : CSR(in_, Qstls(VStoQStlsInput(in_))) { ; }
  
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
