#ifndef DSF_MOMENTS_HPP
#define DSF_MOMENTS_HPP

#include "util/vector2D.hpp"
#include <memory>
#include <vector>

class Input;

namespace thermoUtil {

  /**
   * @brief Compute the first DSF moment using adaptive refinement around
   * under-resolved dielectric peaks.
   */
  std::vector<double> computeAdaptiveFirstMomentFromDSF(
      const std::shared_ptr<const Input> &in,
      const std::vector<double> &wvg,
      const std::vector<double> &frequency,
      const double &mu,
      const Vector2D &lfc,
      const Vector2D &dsf);

  /**
   * @brief Compute the third DSF moment using adaptive refinement around
   * under-resolved dielectric peaks.
   */
  std::vector<double> computeAdaptiveThirdMomentFromDSF(
      const std::shared_ptr<const Input> &in,
      const std::vector<double> &wvg,
      const std::vector<double> &frequency,
      const double &mu,
      const Vector2D &lfc,
      const Vector2D &dsf);

} // namespace thermoUtil

#endif
