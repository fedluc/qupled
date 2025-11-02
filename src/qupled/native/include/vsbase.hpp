#ifndef VSBASE_HPP
#define VSBASE_HPP

#include "input.hpp"
#include "logger.hpp"
#include "mpi_util.hpp"
#include "numerics.hpp"
#include "stls.hpp"
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

  // Constructor
  explicit VSBase() {}
  // Destructor
  virtual ~VSBase() = default;
  // Compute vs-stls scheme
  int compute();
  // Getters
  const std::vector<std::vector<double>> &getFreeEnergyIntegrand() const;
  const std::vector<double> &getFreeEnergyGrid() const;
  double getAlpha() const { return alpha; }

protected:

  // Free parameter
  double alpha;
  // Thermodynamic properties (this must be set from the derived classes)
  std::shared_ptr<ThermoPropBase> thermoProp;
  // Input parameters
  virtual const VSInput &in() const = 0;
  // Compute free parameter
  virtual double computeAlpha() = 0;
  // Initialize
  virtual void init() = 0;
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

class ThermoPropBase : public Logger {

public:

  // Constructor
  explicit ThermoPropBase(const std::shared_ptr<const VSInput> &inPtr_);
  // Destructor
  virtual ~ThermoPropBase() = default;
  // Set the value of the free parameter in the structural properties
  void setAlpha(const double &alpha);
  // Compute the thermodynamic properties
  void compute();
  // Get structural properties
  const std::vector<double> &getSsf();
  const Vector2D &getLfc();
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
  const double &getAlpha() const { return alpha; }

protected:

  using SIdx = StructIdx;
  using Idx = ThermoIdx;
  // Input parameters
  const std::shared_ptr<const VSInput> inPtr;
  // Map between struct and thermo indexes
  static constexpr int NPOINTS = 3;
  // Structural properties (this must be set from the derived classes)
  std::shared_ptr<CSR> structProp;
  // Grid for thermodyamic integration
  std::vector<double> rsGrid;
  // Free energy integrand for NPOINTS state points
  std::vector<std::vector<double>> fxcIntegrand;
  // Free parameter
  double alpha;
  // Flags marking particular state points
  bool isZeroCoupling;
  bool isZeroDegeneracy;
  // Index of the target state point in the free energy integrand
  size_t fxcIdxTargetStatePoint;
  // Access input pointer
  const VSInput &in() const { return *inPtr; }
  // Cast the input member to an Input type
  const Input &inRpa() const {
    return *StlsUtil::dynamic_pointer_cast<VSInput, Input>(inPtr);
  }
  // Compute the free energy
  double computeFreeEnergy(const ThermoPropBase::SIdx iStruct,
                           const bool normalize) const;
  // Build the integration grid
  void setRsGrid();
  // Build the free energy integrand
  void setFxcIntegrand();
  // Set the index of the target state point in the free energy integrand
  void setFxcIdxTargetStatePoint();
  // Get index to acces the structural properties
  ThermoPropBase::SIdx getStructPropIdx();
};

// -----------------------------------------------------------------
// CSR class
// -----------------------------------------------------------------

class CSR {

public:

  // Enumerator to denote the numerical schemes used for the derivatives
  enum Derivative { CENTERED, FORWARD, BACKWARD };
  // Constructor (only manager instances can be created from the outside)
  CSR()
      : CSR(true) {}
  // Destructor
  virtual ~CSR() = default;
  // Solve the scheme
  virtual int compute() = 0;
  // Set the free parameter
  void setAlpha(const double &alpha);
  // Getters
  const std::vector<double> &getSsf(const size_t &idx) const;
  const Vector2D &getLfc(const size_t &idx) const;
  double getCoupling(const size_t &idx) const;
  double getDegeneracy(const size_t &idx) const;
  double getUInt(const size_t &idx) const;
  double getFxcIntegrand(const size_t &idx) const;
  double getAlpha(const size_t &idx) const;
  double getAlpha() { return getAlpha(StructIdx::RS_THETA); }
  double getError() const { return computeError(StructIdx::RS_THETA); }

protected:

  struct DerivativeData {
    Derivative type;
    const Vector2D *up;
    const Vector2D *down;
  };
  // Constructor
  CSR(const bool isManager_)
      : isManager(isManager_),
        isInitialized(false),
        alpha(DEFAULT_ALPHA){};
  // Default value of alpha
  static constexpr double DEFAULT_ALPHA = numUtil::Inf;
  static constexpr int NRS = 3;
  static constexpr int NTHETA = 3;
  // Workers that solve the dielectric scheme
  std::vector<std::shared_ptr<CSR>> workers;
  // Flag marking if this is a manager instance or a worker instance
  const bool isManager;
  // Flag marking if init was already called
  bool isInitialized;
  // Derivative contribution to  the local field correction
  Vector2D lfcDerivative;
  // Free parameter
  double alpha;
  // Data for the local field correction with modified coupling paramter
  DerivativeData lfcRs;
  // Data for the local field correction with modified degeneracy parameter
  DerivativeData lfcTheta;
  // Input data
  virtual const VSInput &inVS() const = 0;
  virtual const Input &inRpa() const = 0;
  // Setup workers
  template <typename Scheme, typename Input>
  void setupWorkers(const Input &in) {
    const double &drs = in.getCouplingResolution();
    const double &dTheta = in.getDegeneracyResolution();
    // If there is a risk of having negative state parameters, shift the
    // parameters so that rs - drs = 0 and/or theta - dtheta = 0
    const double rs = std::max(in.getCoupling(), drs);
    const double theta = std::max(in.getDegeneracy(), dTheta);
    // Setup auxiliary state points
    for (const double &thetaTmp : {theta - dTheta, theta, theta + dTheta}) {
      for (const double &rsTmp : {rs - drs, rs, rs + drs}) {
        std::shared_ptr<Input> inTmp = std::make_shared<Input>(in);
        inTmp->setDegeneracy(thetaTmp);
        inTmp->setCoupling(rsTmp);
        workers.push_back(make_shared<Scheme>(inTmp, false));
      }
    }
    assert(workers.size() == NRS * NTHETA);
    setupDerivativeData();
  }
  void setupDerivativeData();
  // Compute the derivative component of the local field correction
  void computeLfcDerivative();
  // Helper methods to compute the derivatives
  double derivative(const Vector2D &f,
                    const int &l,
                    const size_t &idx,
                    const Derivative &type) const;
  double derivative(const double &f0,
                    const double &f1,
                    const double &f2,
                    const Derivative &type) const;
  // Methods called by compute
  void init();
  void computeLfc();
  void computeSsf();
  void initialGuess();
  void updateSolution();
  double computeError(const size_t &idx) const;
  double computeError() const { return computeError(StructIdx::RS_THETA); }
  virtual void initWorker() = 0;
  virtual void computeLfcWorker() = 0;
  virtual void computeSsfWorker() = 0;
  virtual void initialGuessWorker() = 0;
  virtual void updateSolutionWorker() = 0;
  virtual double computeErrorWorker() const = 0;
  virtual Vector2D &getLfc() = 0;
  // Getters
  virtual const std::vector<double> &getSsf() const = 0;
  virtual const std::vector<double> &getWvg() const = 0;
  virtual const Vector2D &getLfc() const = 0;
  // Convenience methods to handle the workers
  void forEachWorker(auto &&f);
  decltype(auto) withWorker(std::size_t idx, auto &&f) const;
};

#endif
