#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include <limits>
#include <map>
#include "stls.hpp"

// -----------------------------------------------------------------
// Solver for the VS-STLS scheme
// -----------------------------------------------------------------

class StlsCSR : public Stls {
  
  friend class StructProp;
    
private:

  // Default value of alpha
  static constexpr double DEFAULT_ALPHA = numUtil::Inf;
  // Input data
  const VSStlsInput in;
  // Enumerator to denote the numerical schemes used for the derivatives
  enum Derivative { CENTERED,
		    FORWARD,
		    BACKWARD };
  // Stls static local field correction
  vector<double> slfcStls;
  // Free parameter
  double alpha;
  // Pointer to the static local field correction with rs = rs + drs
  vector<double>* slfcStlsRsUp;
  // Pointer to the static local field correction with rs = rs - drs
  vector<double>* slfcStlsRsDown;
  // Pointer to the static local field correction with theta = theta + dtheta
  vector<double>* slfcStlsThetaUp;
  // Pointer to the static local field correction with theta = theta - dtheta
  vector<double>* slfcStlsThetaDown;
  // Numerical scheme used to compute the coupling parameter derivative
  Derivative dTypeRs;
  // Numerical scheme used to compute the degeneracy parametere derivative
  Derivative dTypeTheta;
  // Set the data to compute the coupling parameter derivative
  void setDrsData(vector<double> &slfcStlsRsUp,
		  vector<double> &slfcSltsRsDown,
		  const Derivative &dTypeRs);
  // Set the data to compute the degeneracy parameter derivative
  void setDThetaData(vector<double> &slfcStlsThetaUp,
		     vector<double> &slfcStlsThetaDown,
		     const Derivative &dTypeTheta);
  // Helper methods to compute the derivatives
  double getDerivative(const vector<double>& f,
		       const size_t& idx,
		       const Derivative& type);
  double getDerivative(const double& f0,
		       const double& f1,
		       const double& f2,
		       const Derivative& type);
  // Compute static local field correction
  void computeSlfcStls();
  void computeSlfc();
  // Set the free parameter
  void setAlpha(const double& alpha) { this->alpha = alpha; }
  // Compute the internal energy
  double getInternalEnergy() const;
  // Compute the free energy integrand
  double getFreeEnergyIntegrand() const;
  
public:

  // Constructor
  StlsCSR(const VSStlsInput& in_) : Stls(in_, false, false),
				    in(in_),
				    alpha(DEFAULT_ALPHA),
				    slfcStlsRsUp(nullptr),
				    slfcStlsRsDown(nullptr),
				    slfcStlsThetaUp(nullptr),
				    slfcStlsThetaDown(nullptr),
				    dTypeRs(CENTERED),
				    dTypeTheta(CENTERED) { ; }
  
};


class StructProp {
  
public:

  static constexpr int NPOINTS = 9;
  static constexpr int THETASTEP = 3;
  enum Idx {
    RS_DOWN_THETA_DOWN,
    RS_THETA_DOWN,
    RS_UP_THETA_DOWN,
    RS_DOWN_THETA,
    RS_THETA,
    RS_UP_THETA,
    RS_DOWN_THETA_UP,
    RS_THETA_UP,
    RS_UP_THETA_UP,
  };
  
private:

  // Typdef
  using StlsCSRPtr = std::shared_ptr<StlsCSR>;
  // Vector containing NPOINTS state points to be solved simultaneously
  vector<StlsCSRPtr> stls;
  // Flag marking if the initialization for the stls data was already done
  bool stlsIsInitialized;
  // Perform iterations to compute structural properties
  void doIterations();
  // Generic getter function to return vector data
  const vector<double>& getBase(function<double(const StlsCSRPtr&)> f) const;
  // Vector used as output parameter in the getters functions
  mutable vector<double> outVector;
  
public:

  // Constructor
  StructProp(const VSStlsInput &in_);
  // Compute structural properties
  int compute();
  // Set free parameter
  void setAlpha(const double& alpha);
  // Get coupling parameters for all the state points
  vector<double> getCouplingParameters() const;
  // Get degeneracy parameters for all the state points
  vector<double> getDegeneracyParameters() const;
  // Get internal energy for all the state points
  vector<double> getInternalEnergy() const;
  // Get free energy integrand for all the state points
  vector<double> getFreeEnergyIntegrand() const;
  // Get structural properties for output
  const StlsCSR& getStls(const Idx& idx) const { return *(stls[idx]); }
  // Boolean marking whether the structural properties where computed or not
  bool isComputed() const { return stlsIsInitialized; }
  
};


class ThermoProp {

public:
  
private:

  using SIdx = StructProp::Idx;
  enum Idx {
    THETA_DOWN,
    THETA,
    THETA_UP
  };
  // Map between struct and thermo indexes
  static constexpr int NPOINTS = 3;
  // Structural properties
  StructProp structProp;
  // Grid for thermodyamic integration
  vector<double> rsGrid;
  // Free energy integrand for NPOINTS state points
  vector<vector<double>> fxcIntegrand;
  // Compute the free energy
  double computeFreeEnergy(const SIdx iStruct,
			   const bool normalize) const ;
  
public:

  // Constructors
  ThermoProp(const VSStlsInput &in);
  ThermoProp(const VSStlsInput &in,
	     const ThermoProp &other);
  // Set the value of the free parameter in the structural properties
  void setAlpha(const double& alpha);
  // Compute the thermodynamic properties
  void compute(const VSStlsInput &in);
  const StlsCSR& getStructProp();
  // Get free energy and free energy derivatives
  vector<double> getFreeEnergyData() const;
  // Get internal energy and internal energy derivatives
  vector<double> getInternalEnergyData() const;
  // Get free energy integrand
  const vector<vector<double>>& getFreeEnergyIntegrand() const { return fxcIntegrand; }
  // Get free energy grid
  const vector<double>& getFreeEnergyGrid() const  { return rsGrid; }
  
};

class VSStls : public StlsBase {

private: 

  // Input data
  VSStlsInput in;
  // Thermodynamic properties
  ThermoProp thermoProp;
  // Free parameter
  double alpha;
  // Output verbosity
  const bool verbose;
  // Compute free parameter
  double computeAlpha();
  // Iterations to solve the vs-stls scheme
  void doIterations();
  void updateSolution();
  double alphaDifference(const double& alphaTmp);
 
public:

  // Constructor from initial data
  VSStls(const VSStlsInput &in_) : StlsBase(in_), in(in_),
				   thermoProp(in_),
				   verbose(true) { ; }
  // Constructor for recursive calculations
  VSStls(const VSStlsInput &in_,
	 const ThermoProp& thermoProp_) : StlsBase(in_), in(in_),
					  thermoProp(in_, thermoProp_),
					  verbose(false) { ; }
  // Compute stls scheme
  int compute();
  // Getters
  const ThermoProp& getThermoProp() const { return thermoProp; }
  vector<vector<double>> getFreeEnergyIntegrand() const;
  vector<double> getFreeEnergyGrid() const;
  
};


#endif
