#ifndef THERMO_UTIL_HPP
#define THERMO_UTIL_HPP

#include "util/dimensions_util.hpp"
#include "util/vector2D.hpp"
#include <memory>
#include <vector>

// Forward declarations
class Input;

/**
 * @brief High-level utility functions for computing thermodynamic properties.
 *
 * These free functions provide a convenient interface for evaluating the most
 * commonly requested thermodynamic quantities directly from the arrays returned
 * by the dielectric scheme solvers.
 */
namespace thermoUtil {

  /**
   * @brief Compute the interaction energy per particle.
   *
   * Integrates @f$v(k)\,[S(k) - 1]@f$ over the wave-vector grid using the
   * dimension-appropriate Coulomb potential.
   *
   * @param wvg      Wave-vector grid.
   * @param ssf      Static structure factor values (same length as @p wvg).
   * @param coupling Coupling parameter @f$r_s@f$.
   * @param dim      Spatial dimension.
   * @return Dimensionless interaction energy per particle.
   */
  double computeInternalEnergy(const std::vector<double> &wvg,
                               const std::vector<double> &ssf,
                               const double &coupling,
                               const dimensionsUtil::Dimension &dim);

  /**
   * @brief Compute the exchange-correlation free energy (without
   * normalization).
   *
   * Integrates @f$r_s \cdot u(r_s) / r_s@f$ over the coupling-parameter grid
   * up to @p coupling.
   *
   * @param grid     Coupling-parameter grid.
   * @param rsu      Pre-computed @f$r_s \cdot u(r_s)@f$ values (same length as
   * @p grid).
   * @param coupling Target coupling parameter (upper integration limit).
   * @return Exchange-correlation free energy.
   */
  double computeFreeEnergy(const std::vector<double> &grid,
                           const std::vector<double> &rsu,
                           const double &coupling);

  /**
   * @brief Compute the exchange-correlation free energy with optional
   * normalization.
   *
   * @param grid      Coupling-parameter grid.
   * @param rsu       Pre-computed @f$r_s \cdot u(r_s)@f$ values.
   * @param coupling  Target coupling parameter.
   * @param normalize If true, divide the result by @f$r_s^2@f$.
   * @return Exchange-correlation free energy (optionally normalized).
   */
  double computeFreeEnergy(const std::vector<double> &grid,
                           const std::vector<double> &rsu,
                           const double &coupling,
                           const bool normalize);

  /**
   * @brief Compute the radial distribution function over a real-space grid.
   *
   * Fourier-transforms the static structure factor to obtain g(r) at each
   * point in @p r.
   *
   * @param r    Real-space distance grid.
   * @param wvg  Wave-vector grid.
   * @param ssf  Static structure factor values.
   * @param dim  Spatial dimension.
   * @return Vector of g(r) values, one per entry in @p r.
   */
  std::vector<double> computeRdf(const std::vector<double> &r,
                                 const std::vector<double> &wvg,
                                 const std::vector<double> &ssf,
                                 const dimensionsUtil::Dimension &dim);

  /**
   * @brief Compute the non-interacting (Hartree-Fock) imaginary-time
   * correlation function.
   *
   * Evaluates the HF ITCF for the given wave-vector grid and imaginary-time
   * values using only the ideal density response (no interactions).
   *
   * @param in         Shared pointer to the input parameters.
   * @param wvg        Wave-vector grid.
   * @param tauValues  Imaginary time values in [0, 1] (normalised by beta).
   * @param mu         Chemical potential.
   * @param idr        Ideal density response (rows = wave-vectors, columns =
   *                   Matsubara frequencies).
   * @return 2D vector containing ITCF values (rows = wave-vectors, columns =
   * tau values).
   */
  Vector2D computeItcfNonInteracting(const std::shared_ptr<const Input> &in,
                                     const std::vector<double> &wvg,
                                     const std::vector<double> &tauValues,
                                     const double mu,
                                     const Vector2D &idr);

  /**
   * @brief Compute the imaginary-time correlation function.
   *
   * First computes the non-interacting (HF) contribution. For HF theory
   * (@p in->getTheory() == "HF") the result is returned immediately. For all
   * other theories the interacting correction is applied using the local field
   * correction.
   *
   * @param in         Shared pointer to the input parameters.
   * @param wvg        Wave-vector grid.
   * @param tauValues  Imaginary time values in [0, 1] (normalised by beta).
   * @param mu         Chemical potential.
   * @param idr        Ideal density response (rows = wave-vectors, columns =
   *                   Matsubara frequencies).
   * @param lfc        Local field correction (rows = wave-vectors, columns =
   *                   Matsubara frequencies or 1 for static LFC). Ignored for
   *                   HF theory.
   * @return 2D vector containing ITCF values (rows = wave-vectors, columns =
   * tau values).
   */
  Vector2D computeItcf(const std::shared_ptr<const Input> &in,
                       const std::vector<double> &wvg,
                       const std::vector<double> &tauValues,
                       const double mu,
                       const Vector2D &idr,
                       const Vector2D &lfc);

} // namespace thermoUtil

#endif
