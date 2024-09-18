#ifndef VSBASE_HPP
#define VSBASE_HPP

#include "input.hpp"
#include "logger.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "thermo_util.hpp"
#include "vector2D.hpp"
#include "vector_util.hpp"
#include <limits>
#include <map>
#include <memory>

class ThermoPropBase;
class StructPropBase;
class CSR;

// -----------------------------------------------------------------
// VSBase class
// -----------------------------------------------------------------

class VSBase : public Logger {

public:

  // Constructor from initial data
  explicit VSBase(const VSInput &in_)
      : VSBase(in_, MPIUtil::isRoot()) {}
  VSBase(const VSInput &in_, const bool &verbose_)
      : Logger(verbose_),
        in(in_) {}
  // Destructor
  virtual ~VSBase() = default;
  // Compute vs-stls scheme
  int compute();
  // Getters
  const std::vector<std::vector<double>> &getFreeEnergyIntegrand() const;
  const std::vector<double> &getFreeEnergyGrid() const;
  const std::vector<double> &getAlpha() const;

protected:

  // Input data
  VSInput in;
  // Free parameter
  double alpha;
  // Thermodynamic properties (this must be set from the derived classes)
  std::shared_ptr<ThermoPropBase> thermoProp;
  // Compute free parameter
  virtual double computeAlpha() = 0;
  // Initialize
  void init();
  virtual void initScheme() = 0;
  virtual void initFreeEnergyIntegrand() = 0;
  // Iterations to solve the vs scheme
  void doIterations();
  // Object function used in the secant solver
  double alphaDifference(const double &alphaTmp);
  // Update structural output solution
  virtual void updateSolution() = 0;
};

// -----------------------------------------------------------------
// Indexes
// -----------------------------------------------------------------

enum ThermoIdx { THETA_DOWN, THETA, THETA_UP };

enum StructIdx {
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

// -----------------------------------------------------------------
// ThermoPropBase class
// -----------------------------------------------------------------

class ThermoPropBase {

public:

  // Constructor
  explicit ThermoPropBase(const VSInput &inVS, const RpaInput &inRpa);
  // Destructor
  virtual ~ThermoPropBase() = default;
  // Set the value of the free parameter in the structural properties
  void setAlpha(const double &alpha);
  // Copy free energy integrand
  void copyFreeEnergyIntegrand(const ThermoPropBase &other);
  // Check if there are unsolved state points in the free energy integrand
  bool isFreeEnergyIntegrandIncomplete() const;
  // Get the first unsolved state point in the free energy integrand
  double getFirstUnsolvedStatePoint() const;
  // Compute the thermodynamic properties
  void compute();
  // Get structural properties
  const std::vector<double> &getSsf();
  const std::vector<double> &getSlfc();
  // Get free energy and free energy derivatives
  std::vector<double> getFreeEnergyData() const;
  // Get internal energy and internal energy derivatives
  std::vector<double> getInternalEnergyData() const;
  // Get free energy integrand
  const std::vector<std::vector<double>> &getFreeEnergyIntegrand() const {
    return fxcIntegrand;
  }
  // Get free energy grid
  const std::vector<double> &getFreeEnergyGrid() const { return rsGrid; }
  // Get free parameter values except the last one
  const std::vector<double> &getAlpha() const { return alpha; }

protected:

  using SIdx = StructIdx;
  using Idx = ThermoIdx;
  // Map between struct and thermo indexes
  static constexpr int NPOINTS = 3;
  // Structural properties (this must be set from the derived classes)
  std::shared_ptr<StructPropBase> structProp;
  // Grid for thermodyamic integration
  std::vector<double> rsGrid;
  // Free parameter values for all the coupling parameters stored in rsGrid
  std::vector<double> alpha;
  // Free energy integrand for NPOINTS state points
  std::vector<std::vector<double>> fxcIntegrand;
  // Flags marking particular state points
  bool isZeroCoupling;
  bool isZeroDegeneracy;
  // Index of the target state point in the free energy integrand
  size_t fxcIdxTargetStatePoint;
  // Index of the first unsolved state point in the free energy integrand
  size_t fxcIdxUnsolvedStatePoint;
  // Compute the free energy
  double computeFreeEnergy(const ThermoPropBase::SIdx iStruct,
                           const bool normalize) const;
  // Build the integration grid
  void setRsGrid(const VSInput &inVS, const RpaInput &inRpa);
  // Build the free energy integrand
  void setFxcIntegrand(const VSInput &in);
  // Build the free parameter vector
  void setAlpha(const VSInput &in);
  // Set the index of the target state point in the free energy integrand
  void setFxcIdxTargetStatePoint(const RpaInput &in);
  // Set the index of the first unsolved state point in the free energy
  // integrand
  void setFxcIdxUnsolvedStatePoint();
  // Get index to acces the structural properties
  ThermoPropBase::SIdx getStructPropIdx();
};

// -----------------------------------------------------------------
// StructPropBase class
// -----------------------------------------------------------------

class StructPropBase {

public:

  // Typedef
  using Idx = StructIdx;
  static constexpr int NRS = 3;
  static constexpr int NTHETA = 3;
  static constexpr int NPOINTS = NRS * NTHETA;
  // Constructor
  explicit StructPropBase();
  // Destructor
  virtual ~StructPropBase() = default;
  // Compute structural properties
  int compute();
  // Set free parameter
  void setAlpha(const double &alpha);
  // Get coupling parameters for all the state points
  const std::vector<double> &getCouplingParameters() const;
  // Get degeneracy parameters for all the state points
  const std::vector<double> &getDegeneracyParameters() const;
  // Get internal energy for all the state points
  const std::vector<double> &getInternalEnergy() const;
  // Get free energy integrand for all the state points
  const std::vector<double> &getFreeEnergyIntegrand() const;
  // Get the free parameter
  double getAlpha() const;
  // Get structural properties for output
  const CSR &getCsr(const Idx &idx) const;
  // Boolean marking whether the structural properties where computed or not
  bool isComputed() const { return computed; }

protected:

  // Vector containing NPOINTS state points to be solved simultaneously (this
  // must be defined in the derived class)
  std::vector<std::shared_ptr<CSR>> csr;
  // Flag marking whether the initialization for the stls data is done
  bool csrIsInitialized;
  // Flag marking whether the structural properties were computed
  bool computed;
  // Vector used as output parameter in the getters functions
  mutable std::vector<double> outVector;
  // Setup dependencies for CSR objects
  void setupCSRDependencies();
  // Perform iterations to compute structural properties
  virtual void doIterations() = 0;
  // Generic getter function to return vector data
  const std::vector<double> &
  getBase(std::function<double(const CSR &)> f) const;
};

// -----------------------------------------------------------------
// CSR (compressibility-sum-rule) class
// -----------------------------------------------------------------

class CSR {

public:

  // Enumerator to denote the numerical schemes used for the derivatives
  enum Derivative { CENTERED, FORWARD, BACKWARD };
  // Data for the local field correction with modified state point
  struct DerivativeData {
    Derivative type;
    std::shared_ptr<Vector2D> up;
    std::shared_ptr<Vector2D> down;
  };
  // Constructor
  CSR(const VSInput &inVS_, const RpaInput &inRpa_)
      : inVS(inVS_),
        inRpa(inRpa_),
        lfc(std::make_shared<Vector2D>()),
        alpha(DEFAULT_ALPHA) {}
  // Destructor
  virtual ~CSR() = default;
  // Set the data to compute the coupling parameter derivative
  void setDrsData(CSR &csrRsUp, CSR &csrRsDown, const Derivative &dTypeRs);
  // Set the data to compute the degeneracy parameter derivative
  void setDThetaData(CSR &csrThetaUp,
                     CSR &csrThetaDown,
                     const Derivative &dTypeTheta);
  // Set the free parameter
  void setAlpha(const double &alpha) { this->alpha = alpha; }
  // Get the free parameter
  double getAlpha() const { return alpha; }
  // Get input
  double getCoupling() const { return inRpa.getCoupling(); }
  double getDegeneracy() const { return inRpa.getDegeneracy(); }
  // Compute the internal energy
  double getInternalEnergy() const;
  // Compute the free energy integrand
  double getFreeEnergyIntegrand() const;
  // Publicly esposed private scheme methods
  virtual void init() = 0;
  virtual void initialGuess() = 0;
  virtual void computeSsf() = 0;
  virtual double computeError() = 0;
  virtual void updateSolution() = 0;
  virtual const std::vector<double> &getSsf() const = 0;
  virtual const std::vector<double> &getSlfc() const = 0;
  virtual const std::vector<double> &getWvg() const = 0;

protected:

  // Default value of alpha
  static constexpr double DEFAULT_ALPHA = numUtil::Inf;
  // Input data
  const VSInput inVS;
  const RpaInput inRpa;
  // local field correction (static or dynamic)
  std::shared_ptr<Vector2D> lfc;
  // Free parameter
  double alpha;
  // Data for the local field correction with modified coupling paramter
  DerivativeData lfcRs;
  // Data for the local field correction with modified degeneracy parameter
  DerivativeData lfcTheta;
  // Compute the local field correction with the derivatives contribution
  Vector2D getDerivativeContribution() const;
  // Helper methods to compute the derivatives
  double getDerivative(const std::shared_ptr<Vector2D> &f,
                       const int &l,
                       const size_t &idx,
                       const Derivative &type) const;
  double getDerivative(const double &f0,
                       const double &f1,
                       const double &f2,
                       const Derivative &type) const;
};

#endif
