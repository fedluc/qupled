#ifndef VS_GRID_POINT_HPP
#define VS_GRID_POINT_HPP

#include <cstddef>

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
  inline constexpr GridPoint RS_DOWN_THETA_DOWN = {GridPoint::Rs::DOWN,
                                                   GridPoint::Theta::DOWN};
  inline constexpr GridPoint RS_THETA_DOWN = {GridPoint::Rs::CENTER,
                                              GridPoint::Theta::DOWN};
  inline constexpr GridPoint RS_UP_THETA_DOWN = {GridPoint::Rs::UP,
                                                 GridPoint::Theta::DOWN};
  inline constexpr GridPoint RS_DOWN_THETA = {GridPoint::Rs::DOWN,
                                              GridPoint::Theta::CENTER};
  inline constexpr GridPoint CENTER = {GridPoint::Rs::CENTER,
                                       GridPoint::Theta::CENTER};
  inline constexpr GridPoint RS_UP_THETA = {GridPoint::Rs::UP,
                                            GridPoint::Theta::CENTER};
  inline constexpr GridPoint RS_DOWN_THETA_UP = {GridPoint::Rs::DOWN,
                                                 GridPoint::Theta::UP};
  inline constexpr GridPoint RS_THETA_UP = {GridPoint::Rs::CENTER,
                                            GridPoint::Theta::UP};
  inline constexpr GridPoint RS_UP_THETA_UP = {GridPoint::Rs::UP,
                                               GridPoint::Theta::UP};
} // namespace GridPoints

#endif
