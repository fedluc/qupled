#ifndef INTERNAL_ENERGY_HPP
#define INTERNAL_ENERGY_HPP

#include "util/dimensions_util.hpp"
#include "util/numerics.hpp"
#include <cmath>

/**
 * @brief Computes the interaction (internal) energy per particle.
 *
 * Evaluates the coupling-constant integral
 * @f$u(r_s) = \frac{1}{2} \int_0^\infty v(k)\, [S(k) - 1]\, dk / (2\pi)^d@f$
 * where @f$v(k)@f$ is the Coulomb potential in @f$d@f$ dimensions and
 * @f$S(k)@f$ is the static structure factor. Both 2D and 3D geometries
 * are supported.
 */
class InternalEnergy {

public:

  /**
   * @brief Construct with the integration data and parameters.
   * @param rs_    Coupling parameter (Wigner–Seitz radius).
   * @param yMin_  Lower integration limit (minimum wave-vector).
   * @param yMax_  Upper integration limit (maximum wave-vector).
   * @param ssfi_  Shared pointer to an SSF interpolator.
   * @param itg_   Shared pointer to a 1D numerical integrator.
   * @param dim_   Spatial dimension of the system.
   */
  InternalEnergy(const double &rs_,
                 const double &yMin_,
                 const double &yMax_,
                 std::shared_ptr<Interpolator1D> ssfi_,
                 std::shared_ptr<Integrator1D> itg_,
                 const dimensionsUtil::Dimension &dim_)
      : rs(rs_),
        yMin(yMin_),
        yMax(yMax_),
        itg(itg_),
        ssfi(ssfi_),
        dim(dim_) {}

  /**
   * @brief Compute and return the interaction energy per particle.
   * @return Dimensionless interaction energy @f$u / (e^2 / a_0)@f$.
   */
  double get() const;

private:

  /** @brief Coupling parameter. */
  const double rs;
  /** @brief Lower integration limit (minimum wave-vector). */
  const double yMin;
  /** @brief Upper integration limit (maximum wave-vector). */
  const double yMax;
  /** @brief 1D numerical integrator. */
  const std::shared_ptr<Integrator1D> itg;
  /** @brief Interpolator for the static structure factor. */
  const std::shared_ptr<Interpolator1D> ssfi;
  /** @brief Spatial dimension. */
  const dimensionsUtil::Dimension dim;

  /**
   * @brief Evaluate the integrand at wave-vector @p y.
   * @param y Wave-vector value.
   * @return Integrand value @f$v(y)\,[S(y) - 1]@f$ at @p y.
   */
  double integrand(const double &y) const;

  /**
   * @brief Evaluate the interpolated SSF at wave-vector @p y.
   * @param y Wave-vector value.
   * @return Interpolated SSF value.
   */
  double ssf(const double &y) const;
};

#endif
