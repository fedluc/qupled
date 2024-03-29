#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include <limits>
#include <map>
#include "stls.hpp"

// Forward declarations
class VSStlsInput;

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
  std::vector<double> slfcStls;
  // Free parameter
  double alpha;
  // Pointer to the static local field correction with rs = rs + drs
  std::vector<double>* slfcStlsRsUp;
  // Pointer to the static local field correction with rs = rs - drs
  std::vector<double>* slfcStlsRsDown;
  // Pointer to the static local field correction with theta = theta + dtheta
  std::vector<double>* slfcStlsThetaUp;
  // Pointer to the static local field correction with theta = theta - dtheta
  std::vector<double>* slfcStlsThetaDown;
  // Numerical scheme used to compute the coupling parameter derivative
  Derivative dTypeRs;
  // Numerical scheme used to compute the degeneracy parametere derivative
  Derivative dTypeTheta;
  // Set the data to compute the coupling parameter derivative
  void setDrsData(StlsCSR &stlsRsUp,
		  StlsCSR &sltsRsDown,
		  const Derivative &dTypeRs);
  // Set the data to compute the degeneracy parameter derivative
  void setDThetaData(StlsCSR &stlsThetaUp,
		     StlsCSR &stlsThetaDown,
		     const Derivative &dTypeTheta);
  // Helper methods to compute the derivatives
  double getDerivative(const std::vector<double>& f,
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

  static constexpr int NRS = 3;
  static constexpr int NTHETA = 3;
  static constexpr int NPOINTS = NRS * NTHETA;
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

  // Vector containing NPOINTS state points to be solved simultaneously
  std::vector<StlsCSR> stls;
  // Flag marking whether the initialization for the stls data is done
  bool stlsIsInitialized;
  // Flag marking whether the structural properties were computed
  bool computed;
  // Perform iterations to compute structural properties
  void doIterations();
  // Generic getter function to return vector data
  const std::vector<double>& getBase(std::function<double(const StlsCSR&)> f) const;
  // Vector used as output parameter in the getters functions
  mutable std::vector<double> outVector;
  
public:

  // Constructor
  StructProp(const VSStlsInput &in_);
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
  // Get structural properties for output
  const StlsCSR& getStls(const Idx& idx) const { return stls[idx]; }
  // Boolean marking whether the structural properties where computed or not
  bool isComputed() const { return computed; }
  
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
  std::vector<double> rsGrid;
  // Free energy integrand for NPOINTS state points
  std::vector<std::vector<double>> fxcIntegrand;
  // Flags marking particular state points
  bool isZeroCoupling;
  bool isZeroDegeneracy;
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
  std::vector<double> getFreeEnergyData() const;
  // Get internal energy and internal energy derivatives
  std::vector<double> getInternalEnergyData() const;
  // Get free energy integrand
  const std::vector<std::vector<double>>& getFreeEnergyIntegrand() const { return fxcIntegrand; }
  // Get free energy grid
  const std::vector<double>& getFreeEnergyGrid() const  { return rsGrid; }
  
};

class VSStls : public Rpa {

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
  VSStls(const VSStlsInput &in_) : Rpa(in_, false),
				   in(in_),
				   thermoProp(in_),
				   verbose(true) { ; }
  // Constructor for recursive calculations
  VSStls(const VSStlsInput &in_,
	 const ThermoProp& thermoProp_) : Rpa(in_, false),
					  in(in_),
					  thermoProp(in_, thermoProp_),
					  verbose(false) { ; }
  // Compute vs-stls scheme
  int compute();
  // Getters
  const ThermoProp& getThermoProp() const { return thermoProp; }
  std::vector<std::vector<double>> getFreeEnergyIntegrand() const;
  std::vector<double> getFreeEnergyGrid() const;
  
};


#endif
