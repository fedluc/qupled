#ifndef STATE_POINT_GRID_HPP
#define STATE_POINT_GRID_HPP

#include "input.hpp"
#include "vector2D.hpp"
#include <array>
#include <memory>
#include <vector>

// -----------------------------------------------------------------
// Strong type for addressing a point in the 3x3 state point grid
// -----------------------------------------------------------------

struct GridPoint {

  enum class Rs { DOWN = -1, CENTER = 0, UP = 1 };
  enum class Theta { DOWN = -1, CENTER = 0, UP = 1 };

  Rs rs;
  Theta theta;

  // Maps to flat index 0-8 (theta-outer, rs-inner, matching legacy StructIdx)
  constexpr size_t toIndex() const {
    return static_cast<size_t>((static_cast<int>(theta) + 1) * 3
                               + (static_cast<int>(rs) + 1));
  }

};

// Named constants for all 9 grid points (defined after struct is complete)
namespace GridPoints {
  inline constexpr GridPoint RS_DOWN_THETA_DOWN = {GridPoint::Rs::DOWN,   GridPoint::Theta::DOWN};
  inline constexpr GridPoint RS_THETA_DOWN      = {GridPoint::Rs::CENTER, GridPoint::Theta::DOWN};
  inline constexpr GridPoint RS_UP_THETA_DOWN   = {GridPoint::Rs::UP,     GridPoint::Theta::DOWN};
  inline constexpr GridPoint RS_DOWN_THETA      = {GridPoint::Rs::DOWN,   GridPoint::Theta::CENTER};
  inline constexpr GridPoint CENTER             = {GridPoint::Rs::CENTER, GridPoint::Theta::CENTER};
  inline constexpr GridPoint RS_UP_THETA        = {GridPoint::Rs::UP,     GridPoint::Theta::CENTER};
  inline constexpr GridPoint RS_DOWN_THETA_UP   = {GridPoint::Rs::DOWN,   GridPoint::Theta::UP};
  inline constexpr GridPoint RS_THETA_UP        = {GridPoint::Rs::CENTER, GridPoint::Theta::UP};
  inline constexpr GridPoint RS_UP_THETA_UP     = {GridPoint::Rs::UP,     GridPoint::Theta::UP};
}

// -----------------------------------------------------------------
// StatePointGrid: owns and coordinates 9 perturbed scheme instances
// -----------------------------------------------------------------

template <typename Scheme>
class StatePointGrid {

public:

  using InputType = typename Scheme::InputType;

  explicit StatePointGrid(const std::shared_ptr<const InputType> &in);
  ~StatePointGrid() = default;

  // Set free parameter (propagated to LFC derivative computation)
  void setAlpha(double alpha_) { alpha = alpha_; }
  // Run the full structural property computation for all 9 workers
  int compute();

  // Getters by GridPoint
  const std::vector<double> &getSsf(GridPoint p) const;
  const Vector2D &           getLfc(GridPoint p) const;
  const std::vector<double> &getWvg(GridPoint p) const;
  double getCoupling(GridPoint p)          const;
  double getDegeneracy(GridPoint p)        const;
  double getUInt(GridPoint p)              const;
  double getFxcIntegrandValue(GridPoint p) const;
  double getAlpha() const { return alpha; }
  double getError() const { return lastError; }

  // Direct access to a specific worker (e.g. for QAdder computation)
  const Scheme &getWorkerAt(GridPoint p) const;
  // Direct access to central worker for delegation (idr, sdr, etc.)
  const Scheme &centralWorker() const;

private:

  struct DerivativeData {
    enum class Type { CENTERED, FORWARD, BACKWARD };
    Type   type;
    size_t upIdx;
    size_t downIdx;
  };

  static constexpr int N = 9;

  std::array<std::unique_ptr<Scheme>, N> workers;
  std::array<Vector2D, N>                lfcDerivatives;
  std::array<DerivativeData, N>          rsDerivData;
  std::array<DerivativeData, N>          thetaDerivData;
  std::array<double, N>                  rsValues;
  std::array<double, N>                  thetaValues;
  std::shared_ptr<const InputType>       inPtr;
  double                                 alpha;
  double                                 lastError;
  bool                                   initDone;

  void   setupDerivativeData();
  void   computeLfc();
  void   computeLfcDerivatives();
  double derivative(const Vector2D &f,
                    int l,
                    size_t i,
                    DerivativeData::Type t) const;
  double derivative(double f0,
                    double f1,
                    double f2,
                    DerivativeData::Type t) const;
};

#endif
