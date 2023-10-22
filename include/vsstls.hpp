#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include "stls.hpp"

// -----------------------------------------------------------------
// Solver for the VS-STLS scheme
// -----------------------------------------------------------------

class StlsCSR : public Stls {
  
  friend class StructProp;
    
private:

  // Input data
  const VSStlsInput in;
  // Enumerator to denote the numerical schemes used for the derivatives
  enum Derivative { CENTERED, FORWARD, BACKWARD };
  // Free parameter
  double alpha;
  // Pointer to the static local field correction with rs = rs + drs
  vector<double>* slfcRsUp;
  // Pointer to the static local field correction with rs = rs - drs
  vector<double>* slfcRsDown;
  // Pointer to the static local field correction with theta = theta + dtheta
  vector<double>* slfcThetaUp;
  // Pointer to the static local field correction with theta = theta - dtheta
  vector<double>* slfcThetaDown;
  // Numerical scheme used to compute the coupling parameter derivative
  Derivative dTypeRs;
  // Numerical scheme used to compute the degeneracy parametere derivative
  Derivative dTypeTheta;
  // Set the data to compute the coupling parameter derivative
  void setDrsData(vector<double> &slfcRsUp,
		  vector<double> &slfcRsDown,
		  const Derivative &dTypeRs);
  // Set the data to compute the degeneracy parameter derivative
  void setDThetaData(vector<double> &slfcThetaUp,
		     vector<double> &slfcThetaDown,
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
  void computeSlfc();
  
public:

  // Constructor
  StlsCSR(const VSStlsInput& in_) : Stls(in_, false, false),
				    in(in_),
				    alpha(in_.getAlpha()),
				    slfcRsUp(nullptr),
				    slfcRsDown(nullptr),
				    slfcThetaUp(nullptr),
				    slfcThetaDown(nullptr),
				    dTypeRs(CENTERED),
				    dTypeTheta(CENTERED) { ; }
  // Set derivative data
  void setDerivativeData(std::vector<std::shared_ptr<StlsCSR>>& stlsVector,
			 const size_t& thisIdx);
  // Set the free parameter
  void setAlpha(const double& alpha) { this->alpha = alpha; }
  // Compute the free energy integrand
  double getFreeEnergyIntegrand() const;
  
};


class StructProp {
public:

  static constexpr int NPOINTS = 9;
  static constexpr int STENCIL = 3;
  
private:

  // Input data
  const VSStlsInput in;
  // Vector containing NPOINTS state points to be solved simultaneously
  std::vector<std::shared_ptr<StlsCSR>> stls;
  // Flag marking if the initialization for the stls data was already done
  bool stlsIsInitialized;
  // Perform iterations to compute structural properties
  void doIterations();
  
public:

  // Constructor
  StructProp(const VSStlsInput &in_);
  // Compute structural properties
  int compute();
  // Set free parameter
  void setAlpha(const double& alpha);
  // Get internal energy
  vector<double> getFreeEnergyIntegrand(const double& Theta);
  // Get structural properties for output
  const StlsCSR& getOutputProperties() const { return *(stls[STENCIL + 1]); }
  
};

class VSStls : public StlsBase {

private: 

  // Input data
  VSStlsInput in;
  // Structural properties 
  StructProp structProp;
  // Thermodynamic properties
  std::vector<double> rsGrid;
  vector<vector<double>> freeEnergyIntegrand;
  // ThermoProp internalEnergy;
  // ThermoProp freeEnergy;
  // ThermoProp freeEnergyIntegrand;
  // Free parameter
  double alpha;
  double alphaNew;
  // Private methods
  void init();
  void computeFixedFreeEnergy();
  void doIterations();
  void initialGuess();
  void computeAlpha();
  double computeError();
  void updateSolution();
  void computeFreeEnergyIntegrand();
  
public:

  // Constructors
  VSStls(const VSStlsInput &in_);
  // Compute stls scheme
  int compute();
  
};


#endif
