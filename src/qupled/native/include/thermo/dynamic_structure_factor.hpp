#ifndef DYNAMIC_STRUCTURE_FACTOR_HPP
#define DYNAMIC_STRUCTURE_FACTOR_HPP

#include "util/numerics.hpp"
#include "util/vector2D.hpp"
#include <functional>
#include <memory>
#include <utility>
#include <vector>

class Input;

namespace thermoUtil {

  /**
   * @brief Finite-temperature 3D ideal density response on the real-frequency
   * axis.
   *
   * Computes the real and imaginary parts of Phi_0(x, Omega). The dedicated
   * static expression is used at Omega = 0.
   */
  class IdealDynamicResponse {

  public:

    IdealDynamicResponse(const double &x_,
                         const double &Omega_,
                         const double &Theta_,
                         const double &mu_,
                         std::shared_ptr<Integrator1D> itg_);

    /** @brief Return the real part of Phi_0(x, Omega). */
    double real() const;

    /** @brief Return the imaginary part of Phi_0(x, Omega). */
    double imaginary() const;

  private:

    const double x;
    const double Omega;
    const double Theta;
    const double mu;
    const std::shared_ptr<Integrator1D> itg;

    double realDynamicIntegrand(const double &y) const;
    double realStaticIntegrand(const double &y) const;
    double integrateInfinite(const std::function<double(double)> &func,
                             const std::vector<double> &singularities) const;
  };

  /**
   * @brief Compute the ideal real-frequency response on a rectangular grid.
   *
   * @return Pair containing real and imaginary response arrays, respectively.
   * Rows correspond to wave vectors and columns to frequencies.
   */
  std::pair<Vector2D, Vector2D>
  computeIdealDynamicResponse(const std::vector<double> &wvg,
                              const std::vector<double> &frequency,
                              const double &Theta,
                              const double &mu,
                              const double &intError);

  /**
   * @brief Compute the dynamic structure factor for a static local-field
   * correction using the finite-temperature 3D dielectric response.
   *
   * @param in        Scheme input containing the state point.
   * @param wvg       Positive wave-vector grid.
   * @param frequency Non-negative real-frequency grid.
   * @param mu        Chemical potential.
   * @param lfc       Static local-field correction with shape (nx, 1).
   */
  Vector2D computeDSF(const std::shared_ptr<const Input> &in,
                      const std::vector<double> &wvg,
                      const std::vector<double> &frequency,
                      const double &mu,
                      const Vector2D &lfc);

  /**
   * @brief Evaluate S(x, Omega) for a static local-field correction.
   *
   * This helper is shared by the DSF grid construction and adaptive
   * peak-window quadrature.
   */
  double computeDSFValue(const double &x,
                         const double &Omega,
                         const double &Theta,
                         const double &mu,
                         const double &interaction,
                         const std::shared_ptr<Integrator1D> &itg);

  /**
   * @brief Transform a DSF using adaptive refinement around under-resolved
   * dielectric peaks.
   *
   * Uses CQUAD integration of a cubic DSF interpolant as the baseline, then
   * replaces unresolved peak windows with direct peak-aware quadrature.
   */
  Vector2D computeAdaptiveITCFfromDSF(
      const std::shared_ptr<const Input> &in,
      const std::vector<double> &wvg,
      const std::vector<double> &frequency,
      const std::vector<double> &tau,
      const double &mu,
      const Vector2D &lfc,
      const Vector2D &dsf);

} // namespace thermoUtil

#endif
