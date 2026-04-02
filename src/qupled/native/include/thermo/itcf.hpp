#ifndef ITCF_HPP
#define ITCF_HPP

#include "schemes/input.hpp"
#include "util/dimensions_util.hpp"
#include "util/numerics.hpp"
#include <memory>
#include <span>

namespace thermoUtil {

  /**
   * @brief Computes the non-interacting (Hartree-Fock) imaginary-time
   * correlation function at finite temperature.
   *
   * Evaluates F_HF(x, tau) via numerical integration over the auxiliary
   * momentum. For 3D systems, uses a single 1D integral. For 2D systems,
   * uses a 1D integral with dimension-specific expressions.
   * The SSF is recovered as the special case tau = 0.
   */
  class ItcfNonInteracting : public dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a finite-temperature non-interacting ITCF
     * calculation.
     * @param in_    Shared pointer to the input parameters.
     * @param x_     Wave-vector value.
     * @param mu_    Chemical potential.
     * @param tau_   Imaginary time in [0, 1] (normalised by beta).
     * @param yMin_  Lower integration limit.
     * @param yMax_  Upper integration limit.
     * @param itg_   Shared pointer to a 1D integrator.
     * @param idr0_  Ideal density response at l=0 for the current wave-vector.
     */
    ItcfNonInteracting(const std::shared_ptr<const Input> in_,
                       const double &x_,
                       const double &mu_,
                       const double &tau_,
                       const double &yMin_,
                       const double &yMax_,
                       std::shared_ptr<Integrator1D> itg_,
                       const double &idr0_);

    /**
     * @brief Compute and return the non-interacting ITCF value.
     * @return ITCF value at the current wave-vector and imaginary time.
     */
    double get();

  private:

    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;
    /** @brief Wave-vector. */
    const double x;
    /** @brief Chemical potential. */
    const double mu;
    /** @brief Normalised imaginary time in [0, 1]. */
    const double tau;
    /** @brief Lower integration limit. */
    const double yMin;
    /** @brief Upper integration limit. */
    const double yMax;
    /** @brief 1D numerical integrator. */
    const std::shared_ptr<Integrator1D> itg;
    /** @brief Ideal density response at l=0 for the current wave-vector. */
    const double idr0;
    /** @brief Stores the ITCF result. */
    double res;

    /**
     * @brief Compute the ITCF for 3D systems.
     *
     * Evaluates the ITCF via 1D integration over the auxiliary momentum.
     */
    void compute3D() override;
    /**
     * @brief Compute the ITCF for 2D systems.
     *
     * Evaluates the ITCF via 1D integration over the auxiliary momentum.
     */
    void compute2D() override;
    /**
     * @brief 3D integrand over auxiliary momentum @p y.
     *
     * Evaluates the ITCF integrand for 3D systems at the specified tau.
     * @param y Auxiliary momentum variable.
     */
    double integrand(const double &y) const;
    /**
     * @brief 2D integrand over auxiliary momentum @p y.
     *
     * Evaluates the ITCF integrand for 2D systems at the specified tau.
     * @param y Outer integration variable.
     */
    double integrand2D(const double &y) const;
  };

  /**
   * @brief Helper class for computing the normalized interaction potential.
   */
  class ItcfHelper {

  public:

    /**
     * @brief Construct with input and wave-vector.
     * @param in_ Shared pointer to the input parameters.
     * @param x_  Wave-vector value.
     */
    ItcfHelper(const std::shared_ptr<const Input> in_, const double &x_)
        : in(in_),
          x(x_) {}

    /** @brief Normalized interaction potential at the current wave-vector. */
    double ip() const;

  private:

    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;
    /** @brief Wave-vector value. */
    const double x;
  };

  /**
   * @brief Computes the RPA imaginary-time correlation function at finite
   * temperature.
   *
   * Evaluates F(x, tau) = F_HF(x, tau) minus the Matsubara correction sum
   * weighted by cos(2*pi*l*tau). The SSF is recovered as the special case
   * tau = 0.
   */
  class Itcf : public dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a finite-temperature RPA ITCF calculation.
     * @param x_      Wave-vector value.
     * @param itcfHF_ HF imaginary-time correlation function at this wave-vector
     *                and imaginary time tau.
     * @param lfc_    Span over the local field correction array.
     * @param in_     Shared pointer to the input parameters.
     * @param idr_    Span over the ideal density response array.
     * @param tau_    Imaginary time in [0, 1] (normalised by beta).
     */
    Itcf(const double &x_,
         const double &itcfHF_,
         std::span<const double> lfc_,
         const std::shared_ptr<const Input> in_,
         std::span<const double> idr_,
         const double &tau_);

    /** @brief Compute and return the RPA ITCF value. */
    double get();

  private:

    /** @brief Wave-vector value. */
    const double x;
    /** @brief Hartree-Fock contribution. */
    const double itcfHF;
    /** @brief Local field correction values. */
    std::span<const double> lfc;
    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;
    /** @brief Ideal density response values over Matsubara frequencies. */
    const std::span<const double> idr;
    /** @brief Normalised imaginary time in [0, 1]. */
    const double tau;
    /** @brief Stores the result of the Matsubara summation. */
    double res;

    /**
     * @brief Compute the ITCF for 3D systems.
     *
     * Evaluates F(x, tau) = F_HF(x, tau) - 1.5 * v(x) * Theta * sum_l,
     * where sum_l is the Matsubara frequency summation.
     */
    void compute3D() override;
    /**
     * @brief Compute the ITCF for 2D systems.
     *
     * Evaluates F(x, tau) = F_HF(x, tau) - v(x) * Theta * sum_l,
     * where sum_l is the Matsubara frequency summation.
     */
    void compute2D() override;
    /**
     * @brief Compute the Matsubara frequency summation.
     * @return Sum over Matsubara frequencies weighted by cos(2*pi*l*tau).
     */
    double computeMatsubaraSummation() const;
  };

} // namespace thermoUtil

#endif
