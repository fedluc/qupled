#ifndef DSF_INTEGRATION_HPP
#define DSF_INTEGRATION_HPP

#include "util/vector2D.hpp"
#include <functional>
#include <memory>
#include <vector>

class Input;

namespace thermoUtil {

  /**
   * @brief Numerical controls for adaptive positive-frequency DSF integrals.
   */
  struct DSFIntegrationOptions {
    double peakSearchXMax = 3.0;
    // Largest x=q/k_F where collective-mode peak corrections are attempted.

    double peakSearchOmegaBuffer = 5.0;
    // Extra frequency window above the recoil scale x^2 used when scanning
    // for dielectric roots: Omega_scan <= min(Omega_max, x^2 + buffer).

    double peakSearchStepMax = 0.1;
    // Maximum frequency step used during the dielectric-root sign-change scan.

    int peakSearchMinIntervals = 200;
    // Minimum number of scan intervals used between 0 and the scan cutoff.

    double slopeStepAbs = 1.0e-7;
    // Absolute lower bound for the finite-difference step used to estimate
    // dD_R/dOmega at a collective-mode root.

    double slopeStepRel = 1.0e-5;
    // Relative finite-difference step, scaled by the collective-mode frequency.

    double peakWindowSpacingFactor = 8.0;
    // Minimum peak-replacement half-window measured in frequency-grid spacings.

    double peakWindowWidthFactor = 32.0;
    // Minimum peak-replacement half-window measured in estimated peak widths.

    double narrowPeakSpacingFactor = 1.0e-4;
    // Use the analytic narrow-peak area when the width is this small compared
    // with the frequency-grid spacing.

    double narrowPeakOmegaFactor = 1.0e-12;
    // Use the analytic narrow-peak area when the width is this small compared
    // with the peak frequency scale.
  };

  /**
   * @brief Compute weighted positive-frequency DSF integrals with adaptive
   * refinement around under-resolved dielectric peaks.
   *
   * The stored DSF grid is integrated with CQUAD as a baseline. Collective
   * mode windows are replaced by direct dielectric-response quadrature when
   * the mode width is not sufficiently resolved by the frequency grid.
   */
  Vector2D computeAdaptiveWeightedDSFIntegral(
      const std::shared_ptr<const Input> &in,
      const std::vector<double> &wvg,
      const std::vector<double> &frequency,
      const double &mu,
      const Vector2D &lfc,
      const Vector2D &dsf,
      const std::vector<std::function<double(double)>> &kernels,
      const double &maxWidthFactor,
      const DSFIntegrationOptions &options = DSFIntegrationOptions());

} // namespace thermoUtil

#endif
