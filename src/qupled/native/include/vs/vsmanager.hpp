#ifndef VS_VSMANAGER_HPP
#define VS_VSMANAGER_HPP

#include "dimensions_util.hpp"
#include "vs/grid_point.hpp"
#include "vs/vsworker.hpp"
#include <array>
#include <memory>
#include <vector>

// -----------------------------------------------------------------
// VSManager: manages 9 workers and derivative bookkeeping
// -----------------------------------------------------------------

class VSManager {
public:

  void setAlpha(double alpha_) { alpha = alpha_; }
  double getAlpha() const { return alpha; }
  double getError() const { return lastError; }
  // Virtual compute method (implemented by derived classes)
  virtual int compute() = 0;
  // Getters by GridPoint
  const std::vector<double> &getSsf(GridPoint p) const;
  const Vector2D &getLfc(GridPoint p) const;
  const std::vector<double> &getWvg(GridPoint p) const;
  double getCoupling(GridPoint p) const;
  double getDegeneracy(GridPoint p) const;
  double getUInt(GridPoint p) const;
  double getQAdder(GridPoint p) const;
  double getFxcIntegrandValue(GridPoint p) const;
  const VSWorker &getWorkerAt(GridPoint p) const;
  // Convenience getters for central worker (delegate to worker interface)
  const std::vector<double> &getSsf() const;
  const Vector2D &getLfc() const;
  const std::vector<double> &getWvg() const;
  const Vector2D &getIdr() const;
  std::vector<double> getSdr() const;
  double getUInt() const;

protected:

  VSManager(double drs_,
            double dTheta_,
            double dx_,
            dimensionsUtil::Dimension dim_);

  struct DerivativeData {
    enum class Type { CENTERED, FORWARD, BACKWARD };
    Type type;
    size_t upIdx;
    size_t downIdx;
  };

  static constexpr int N = 9;
  std::array<std::unique_ptr<VSWorker>, N> workers;
  std::array<Vector2D, N> lfcDerivatives;
  std::array<DerivativeData, N> rsDerivData;
  std::array<DerivativeData, N> thetaDerivData;
  std::array<double, N> rsValues;
  std::array<double, N> thetaValues;
  double alpha;
  mutable double lastError;
  bool initDone;
  double drs;
  double dTheta;
  double dx;
  dimensionsUtil::Dimension dim;

  void setupDerivativeData();

  // Non-virtual iteration helpers — subclasses delegate to these
  void init();
  void computeLfc();
  void computeSsf();
  double computeError() const;
  void updateSolution();
  void initialGuess();

private:

  void computeLfcDerivatives();
  double
  derivative(const Vector2D &f, int l, size_t i, DerivativeData::Type t) const;
  double
  derivative(double f0, double f1, double f2, DerivativeData::Type t) const;
};

#endif
