#ifndef VS_VS_MASTER_BASE_HPP
#define VS_VS_MASTER_BASE_HPP

#include "dimensions_util.hpp"
#include "vs/grid_point.hpp"
#include "vs/vsworker_base.hpp"
#include <array>
#include <memory>
#include <vector>

// -----------------------------------------------------------------
// VSMasterBase: manages 9 workers and derivative bookkeeping
// -----------------------------------------------------------------

class VSMasterBase {
public:

  void setAlpha(double alpha_) { alpha = alpha_; }
  double getAlpha() const { return alpha; }
  double getError() const { return lastError; }
  // Getters by GridPoint
  const std::vector<double> &getSsf(GridPoint p) const;
  const Vector2D &getLfc(GridPoint p) const;
  const std::vector<double> &getWvg(GridPoint p) const;
  double getCoupling(GridPoint p) const;
  double getDegeneracy(GridPoint p) const;
  double getUInt(GridPoint p) const;
  double getFxcIntegrandValue(GridPoint p) const;
  const VSWorkerBase &getWorkerAt(GridPoint p) const;

protected:

  VSMasterBase(double drs_,
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
  std::array<std::unique_ptr<VSWorkerBase>, N> workers;
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
  void computeSynchronizedLfc();

  // Non-virtual iteration helpers — subclasses delegate to these
  void masterInit();
  void masterComputeLfc();
  void masterComputeSsf();
  double masterComputeError() const;
  void masterUpdateSolution();
  void masterInitialGuess();

private:

  void computeLfcDerivatives();
  double
  derivative(const Vector2D &f, int l, size_t i, DerivativeData::Type t) const;
  double
  derivative(double f0, double f1, double f2, DerivativeData::Type t) const;
};

#endif
