#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include "stls.hpp"

// -----------------------------------------------------------------
// Solver for the VS-STLS scheme
// -----------------------------------------------------------------

class StlsCSR : public Stls {
  
public:

  // Enumerator to denote the numerical schemes used for the derivatives
  enum Derivative { CENTERED, FORWARD, BACKWARD };
  // Enumerator to denote the actions to perform when compute is called
  enum Action { INITIALIZE, GUESS, SOLUTION, ERROR, UPDATE };
    
private:

  // Input data
  const VSStlsInput in;
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
  void setDerivativeData(std::vector<std::unique_ptr<StlsCSR>>& stlsVector,
			 const size_t& thisIdx);
  // Set the free parameter
  void setAlpha(const double& alpha) { this->alpha = alpha; }
  // Perform the action specified in input
  double doAction(const Action& action);
  
};


class StructProp {
public:

  static constexpr int NPOINTS = 9;
  static constexpr int STENCIL = 3;
  
private:

  // Input data
  const VSStlsInput in;
  // Vector containing NPOINTS state points to be solved simultaneously
  std::vector<std::unique_ptr<StlsCSR>> stls;
  // Perform iterations to compute structural properties
  void doIterations();
  
public:

  // Constructor
  StructProp(const VSStlsInput &in_);
  // Compute structural properties
  int compute();
  // Set free parameter
  void setAlpha(const double& alpha);
  
};

class ThermoProp {
private:

  // Grid used for thermodynamic integration
  std::vector<double> grid;
  // Thermodynamic properties at NPOINTS state points
  std::vector<std::vector<double>> prop;
  
public:
  
  static constexpr int NPOINTS = 3;
  ThermoProp(const VSStlsInput &in);
  
};

class VSStls {

private: 


  // Private members
  VSStlsInput in;
  StructProp structProp;
  ThermoProp freeEnergyIntegrand;
  ThermoProp rsa;
  std::vector<double> freeEnergy;
  std::vector<double> rsCutOff;
  double alpha;
  double alphaNew;
  double aMix;
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
  VSStls(const VSStlsInput &in_) : in(in_), structProp(in_),
				   freeEnergyIntegrand(in_),
				   rsa(in_),
				   freeEnergy(StructProp::NPOINTS),
				   rsCutOff(StructProp::NPOINTS) { ; }
  // Compute stls scheme
  int compute();
  
};


#endif
