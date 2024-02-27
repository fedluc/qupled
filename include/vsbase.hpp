#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include <limits>
#include <map>

// Forward declarations
class VSStlsInput;

// -----------------------------------------------------------------
// Solver for the VS-STLS scheme
// -----------------------------------------------------------------

template <typename T, Scheme>
class CSR : public Scheme {

protected:
  
  // Default value of alpha
  static constexpr double DEFAULT_ALPHA = numUtil::Inf;
  // Enumerator to denote the numerical schemes used for the derivatives
  enum Derivative { CENTERED,
		    FORWARD,
		    BACKWARD };
  // Input data
  const VSStlsInput in;
  // local field correction (static or dynamic)
  T lfc;
  // Free parameter
  double alpha;
  // Pointer to the local field correction with rs = rs + drs
  T* lfcRsUp;
  // Pointer to the local field correction with rs = rs - drs
  T* lfcRsDown;
  // Pointer to the local field correction with theta = theta + dtheta
  T* lfcThetaUp;
  // Pointer to the local field correction with theta = theta - dtheta
  T* lfcThetaDown;
  // Numerical scheme used to compute the coupling parameter derivative
  Derivative dTypeRs;
  // Numerical scheme used to compute the degeneracy parametere derivative
  Derivative dTypeTheta;
  // Set the data to compute the coupling parameter derivative
  void setDrsData(CSR<T> &csrRsUp,
		  CSR<T> &csrRsDown,
		  const Derivative &dTypeRs);
  // Set the data to compute the degeneracy parameter derivative
  void setDThetaData(CSR<T> &csrThetaUp,
		     CSR<T> &csrThetaDown,
		     const Derivative &dTypeTheta);
  // Helper methods to compute the derivatives
  double getDerivative(const T& f,
		       const size_t& idx,
		       const Derivative& type);
  double getDerivative(const double& f0,
		       const double& f1,
		       const double& f2,
		       const Derivative& type);
  // Set the free parameter
  void setAlpha(const double& alpha) { this->alpha = alpha; }
  // Compute the internal energy
  double getInternalEnergy() const;
  // Compute the free energy integrand
  double getFreeEnergyIntegrand() const;
  
public:

  // Constructor
  CSR(const VSStlsInput& in_,
      const bool&&... args) : Scheme(in_, args),
			      in(in_),
			      alpha(DEFAULT_ALPHA),
			      lfcRsUp(nullptr),
			      lfcRsDown(nullptr),
			      lfcThetaUp(nullptr),
			      lfcThetaDown(nullptr),
			      dTypeRs(CENTERED),
			      dTypeTheta(CENTERED) { ; }
};

class StlsCSR : public CSR<std::vector<double>, Stls> {
  
  friend class StructProp;
    
private:

  using CSR = CSR<std::vector<double>, Stls>;
  // Compute static local field correction
  void computeSlfcStls();
  void computeSlfc();
  // Compute the internal energy
  double getInternalEnergy() const;
  // Compute the free energy integrand
  double getFreeEnergyIntegrand() const;
  
public:

  // Constructor
  StlsCSR(const VSStlsInput& in_) : CSR(in_, false, false, int) { ; }
  
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

protected:
  
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

private:
  
  // Structural properties
  StructProp structProp;
  
protected:

  using SIdx = StructProp::Idx;
  enum Idx {
    THETA_DOWN,
    THETA,
    THETA_UP
  };
  // Map between struct and thermo indexes
  static constexpr int NPOINTS = 3;
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

  // Thermodynamic properties
  ThermoProp thermoProp;
  
protected: 

  // Input data
  VSStlsInput in;
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
				   thermoProp(in_),
				   in(in_),
				   verbose(true) { ; }
  // Constructor for recursive calculations
  VSStls(const VSStlsInput &in_,
	 const ThermoProp& thermoProp_) : Rpa(in_, false),
					  thermoProp(in_, thermoProp_),
					  in(in_),
					  verbose(false) { ; }
  // Compute vs-stls scheme
  int compute();
  // Getters
  const ThermoProp& getThermoProp() const { return thermoProp; }
  std::vector<std::vector<double>> getFreeEnergyIntegrand() const;
  std::vector<double> getFreeEnergyGrid() const;
  
};


#endif
