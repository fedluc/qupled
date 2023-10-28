#ifndef VSSTLS_HPP
#define VSSTLS_HPP

#include <limits>
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
  
public:

  // Constructor
  StlsCSR(const VSStlsInput& in_) : Stls(in_, false, false),
				    in(in_),
				    alpha(in_.getAlpha()),
				    slfcStlsRsUp(nullptr),
				    slfcStlsRsDown(nullptr),
				    slfcStlsThetaUp(nullptr),
				    slfcStlsThetaDown(nullptr),
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
  vector<double> getFreeEnergyIntegrand(const double& theta);
  // Get structural properties for output
  const StlsCSR& getStls(const double& rs,
			 const double& theta) const;
  bool isComputed() const { return stlsIsInitialized; }
  
};

class VSStls : public StlsBase {

private: 

  // Typedef
  double Inf = numeric_limits<double>::infinity();  
  using doubleVector = std::vector<std::vector<double>>;
  // Input data
  VSStlsInput in;
  // Structural properties 
  StructProp structProp;
  // Thermodynamic properties
  std::shared_ptr<std::vector<double>> rsGrid;
  std::shared_ptr<doubleVector> fxcIntegrand;
  std::vector<double> rsGridLocal;
  std::vector<double> fxcIntegrandLocal;
  // Free parameter
  double alpha;
  double alphaNew;
  // Output verbosity
  const bool verbose;
  // Initialize data 
  void init();
  void setFreeEnergyIntegrand();
  // Compute free parameter
  void computeAlpha();
  // Compute free energy
  void computeFreeEnergyIntegrand();
  double computeFreeEnergy(const double& rs,
			   const bool& normalize);
  // Iterations to solve the vs-stls scheme
  void doIterations();
  void initialGuess();
  double computeError();
  void updateSolution();
  // Private constructor for nested calculations
  VSStls(const VSStlsInput &in_,
	 std::shared_ptr<std::vector<double>> &rsGrid_,
	 std::shared_ptr<doubleVector> &fxcIntegrand_);
  
public:

  // Constructors
  VSStls(const VSStlsInput &in_);
  // Compute stls scheme
  int compute();
  // Getters
  std::vector<double> getFreeEnergyIntegrand() const {
    return fxcIntegrand->at(0);
  }
  std::vector<double> getFreeEnergyGrid() const { return *rsGrid; }
  
};


#endif
