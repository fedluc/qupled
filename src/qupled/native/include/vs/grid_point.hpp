#ifndef VS_GRID_POINT_HPP
#define VS_GRID_POINT_HPP

#include <cstddef>

/**
 * @brief Strong type for addressing a point in the 3×3 (rs, Theta) state-point
 * grid.
 *
 * The variational-swarm optimization evaluates the free energy on a 3×3 grid
 * of (rs, Theta) values centred on the target state point. @p GridPoint
 * encodes a grid position using typed enumerators for the coupling-parameter
 * (@p Rs) and degeneracy-parameter (@p Theta) axes, and provides a mapping to
 * the flat 0–8 index used internally by @p VSManager.
 */
struct GridPoint {

  /** @brief Position along the coupling-parameter (rs) axis. */
  enum class Rs { DOWN = -1, CENTER = 0, UP = 1 };

  /** @brief Position along the degeneracy-parameter (Theta) axis. */
  enum class Theta { DOWN = -1, CENTER = 0, UP = 1 };

  /** @brief Rs-axis position of this grid point. */
  Rs rs;
  /** @brief Theta-axis position of this grid point. */
  Theta theta;

  /**
   * @brief Map this grid point to a flat index in [0, 8].
   *
   * The mapping uses theta as the outer (slow) index and rs as the inner
   * (fast) index, matching the legacy @p StructIdx convention.
   *
   * @return Flat index in [0, 8].
   */
  constexpr size_t toIndex() const {
    return static_cast<size_t>((static_cast<int>(theta) + 1) * 3
                               + (static_cast<int>(rs) + 1));
  }
};

/**
 * @brief Named constants for all nine grid points in the 3×3 (rs, Theta) grid.
 *
 * The naming convention is RS_{rs-position}_THETA_{theta-position}, where
 * positions are DOWN, (CENTER omitted), or UP.  The centre point is simply
 * @p CENTER.
 */
namespace GridPoints {
  /** @brief Grid point at (rs–, Theta–). */
  inline constexpr GridPoint RS_DOWN_THETA_DOWN = {GridPoint::Rs::DOWN,
                                                   GridPoint::Theta::DOWN};
  /** @brief Grid point at (rs, Theta–). */
  inline constexpr GridPoint RS_THETA_DOWN = {GridPoint::Rs::CENTER,
                                              GridPoint::Theta::DOWN};
  /** @brief Grid point at (rs+, Theta–). */
  inline constexpr GridPoint RS_UP_THETA_DOWN = {GridPoint::Rs::UP,
                                                 GridPoint::Theta::DOWN};
  /** @brief Grid point at (rs–, Theta). */
  inline constexpr GridPoint RS_DOWN_THETA = {GridPoint::Rs::DOWN,
                                              GridPoint::Theta::CENTER};
  /** @brief Grid point at the target state point (rs, Theta). */
  inline constexpr GridPoint CENTER = {GridPoint::Rs::CENTER,
                                       GridPoint::Theta::CENTER};
  /** @brief Grid point at (rs+, Theta). */
  inline constexpr GridPoint RS_UP_THETA = {GridPoint::Rs::UP,
                                            GridPoint::Theta::CENTER};
  /** @brief Grid point at (rs–, Theta+). */
  inline constexpr GridPoint RS_DOWN_THETA_UP = {GridPoint::Rs::DOWN,
                                                 GridPoint::Theta::UP};
  /** @brief Grid point at (rs, Theta+). */
  inline constexpr GridPoint RS_THETA_UP = {GridPoint::Rs::CENTER,
                                            GridPoint::Theta::UP};
  /** @brief Grid point at (rs+, Theta+). */
  inline constexpr GridPoint RS_UP_THETA_UP = {GridPoint::Rs::UP,
                                               GridPoint::Theta::UP};
} // namespace GridPoints

#endif
