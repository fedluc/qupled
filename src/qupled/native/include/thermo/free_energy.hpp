#ifndef FREE_ENERGY_HPP
#define FREE_ENERGY_HPP

#include "util/numerics.hpp"

/**
 * @brief Computes the exchange-correlation free energy by integrating over the
 * coupling parameter.
 *
 * The free energy is obtained from the coupling-constant integration
 * @f$f_{xc}(r_s) = \int_0^{r_s} \frac{u(r_s')}{r_s'}\, dr_s'@f$
 * where @f$u(r_s')@f$ is the interaction energy per particle at coupling
 * parameter @f$r_s'@f$. The integrand @f$r_s' u(r_s')@f$ is provided via a
 * 1D spline interpolator.
 */
class FreeEnergy {

public:

  /**
   * @brief Construct with the integration data and options.
   * @param rs_        Target coupling parameter (upper integration limit).
   * @param rsui_      Interpolator for the integrand @f$r_s \cdot u(r_s)@f$.
   * @param itg_       Shared pointer to a 1D numerical integrator.
   * @param normalize_ If true, divide the result by @f$r_s^2@f$.
   */
  FreeEnergy(const double &rs_,
             std::shared_ptr<Interpolator1D> rsui_,
             std::shared_ptr<Integrator1D> itg_,
             const bool normalize_)
      : rs(rs_),
        itg(itg_),
        rsui(rsui_),
        normalize(normalize_) {}

  /**
   * @brief Compute and return the free energy.
   * @return Exchange-correlation free energy (optionally normalized by
   * @f$r_s^2@f$).
   */
  double get() const;

private:

  /** @brief Target coupling parameter (upper integration limit). */
  const double rs;

  /** @brief 1D numerical integrator. */
  const std::shared_ptr<Integrator1D> itg;

  /** @brief Interpolator for the integrand @f$r_s \cdot u(r_s)@f$. */
  const std::shared_ptr<Interpolator1D> rsui;

  /**
   * @brief Evaluate the coupling-constant integrand at @p y.
   * @param y Coupling parameter value.
   * @return Integrand value @f$(r_s' \cdot u(r_s')) / r_s'@f$ at @p y.
   */
  double integrand(const double y) const;

  /** @brief If true, the result is divided by @f$r_s^2@f$. */
  const bool normalize;
};

#endif
