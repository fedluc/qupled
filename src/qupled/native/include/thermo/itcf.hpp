#ifndef ITCF_HPP
#define ITCF_HPP

#include "schemes/input.hpp"
#include "util/dimensions_util.hpp"
#include "util/numerics.hpp"
#include <memory>
#include <span>

namespace thermoUtil {

  /**
   * @brief Base class holding shared state for ITCF calculations.
   *
   * Provides common state (wave-vector, input, imaginary time) and the
   * normalized interaction potential for both non-interacting and interacting
   * ITCF computations.
   */
  class ItcfBase {

  protected:

    /**
     * @brief Construct the base with quantities needed for ITCF evaluation.
     * @param x_   Wave-vector value.
     * @param in_  Shared pointer to the input parameters.
     * @param tau_ imaginary time.
     */
    ItcfBase(const double &x_,
             const std::shared_ptr<const Input> in_,
             const double &tau_)
        : x(x_),
          in(in_),
          tau(tau_) {}

    /** @brief Wave-vector value. */
    const double x;

    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;

    /** @brief Imaginary time. */
    const double tau;

    /** @brief Normalized interaction potential at the current wave-vector. */
    double ip() const;
  };

  /**
   * @brief Computes the non-interacting (Hartree-Fock) imaginary-time
   * correlation function at finite temperature.
   *
   * Evaluates F_HF(x, tau) via numerical integration over the auxiliary
   * momentum. For 3D systems, uses a single 1D integral. For 2D systems,
   * uses a 1D integral with dimension-specific expressions.
   * The SSF is recovered as the special case tau = 0.
   */
  class ItcfNonInteracting : public ItcfBase,
                             public dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a finite-temperature non-interacting ITCF
     * calculation.
     * @param x_     Wave-vector value.
     * @param in_    Shared pointer to the input parameters.
     * @param tau_   Imaginary time.
     * @param mu_    Chemical potential.
     * @param yMin_  Lower integration limit.
     * @param yMax_  Upper integration limit.
     * @param itg_   Shared pointer to a 1D integrator.
     * @param idr0_  Ideal density response at l=0 for the current wave-vector.
     */
    ItcfNonInteracting(const double &x_,
                       const std::shared_ptr<const Input> in_,
                       const double &tau_,
                       const double &mu_,
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

    /** @brief Chemical potential. */
    const double mu;
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
    double integrand3D(const double &y) const;
    /**
     * @brief 2D integrand over auxiliary momentum @p y.
     *
     * Evaluates the ITCF integrand for 2D systems at the specified tau.
     * @param y Outer integration variable.
     */
    double integrand2D(const double &y) const;
  };

  /**
   * @brief Computes the non-interacting (Hartree-Fock) imaginary-time
   * correlation function at zero temperature.
   *
   * Evaluates F_HF(x, tau) via an analytic formula in 3D and a 2D
   * k-space integral in 2D with constant limits y in [0, 1], phi in [0, pi],
   * and a Heaviside factor in the angular integrand.
   * The SSF is recovered as the special case tau = 0.
   */
  class ItcfNonInteractingGround : public ItcfBase,
                                   public dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a ground-state non-interacting ITCF calculation.
     * @param x_   Wave-vector value.
     * @param in_  Shared pointer to the input parameters.
     * @param tau_ Imaginary time.
     * @param itg2_ Shared pointer to a 2D integrator used by the 2D branch.
     */
    ItcfNonInteractingGround(const double &x_,
                             const std::shared_ptr<const Input> in_,
                             const double &tau_,
                             std::shared_ptr<Integrator2D> itg2_);

    /**
     * @brief Compute and return the non-interacting ground-state ITCF value.
     * @return ITCF value at the current wave-vector and imaginary time.
     */
    double get();

  private:

    /** @brief 2D numerical integrator for the 2D ground-state branch. */
    const std::shared_ptr<Integrator2D> itg2;

    /** @brief Stores the ITCF result. */
    double res;

    /** @brief Compute the ITCF for 3D systems. */
    void compute3D() override;

    /** @brief Compute the ITCF for 2D systems. */
    void compute2D() override;

    /**
     * @brief 2D inner angular integrand for ground-state ITCF.
     *
     * Returns zero when the particle-hole phase-space condition is not met,
     * i.e. when y^2 + x^2 + 2xy cos(phi) - 1 < 0.
     * @param p Angular variable phi.
     */
    double integrand2DIn(const double &p) const;
  };

  /**
   * @brief Computes the RPA imaginary-time correlation function at finite
   * temperature.
   *
   * Evaluates F(x, tau) = F_HF(x, tau) minus the Matsubara correction sum
   * weighted by cos(2*pi*l*tau*Theta), where tau is the imaginary time
   * and Theta is the degeneracy parameter. The SSF is recovered as the special
   * case tau = 0.
   */
  class Itcf : public ItcfBase, public dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a finite-temperature RPA ITCF calculation.
     * @param x_      Wave-vector value.
     * @param in_     Shared pointer to the input parameters.
     * @param tau_    Imaginary time.
     * @param itcfHF_ HF imaginary-time correlation function at this wave-vector
     *                and imaginary time tau.
     * @param lfc_    Span over the local field correction array.
     * @param idr_    Span over the ideal density response array.
     */
    Itcf(const double &x_,
         const std::shared_ptr<const Input> in_,
         const double &tau_,
         const double &itcfHF_,
         std::span<const double> lfc_,
         std::span<const double> idr_);

    /** @brief Compute and return the RPA ITCF value. */
    double get();

  private:

    /** @brief Hartree-Fock contribution. */
    const double itcfHF;
    /** @brief Local field correction values. */
    std::span<const double> lfc;
    /** @brief Ideal density response values over Matsubara frequencies. */
    const std::span<const double> idr;
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

  /**
   * @brief Computes the RPA imaginary-time correlation function at zero
   * temperature.
   *
   * Evaluates F(x, tau) = F_HF(x, tau) + (3/2pi) * integral over real
   * frequencies of the RPA correction weighted by exp(-Omega*tau).
   * The SSF is recovered as the special case tau = 0.
   */
  class ItcfGround : public ItcfBase {

  public:

    /**
     * @brief Construct for a ground-state RPA ITCF calculation.
     * @param x_      Wave-vector value.
     * @param in_     Shared pointer to the input parameters.
     * @param tau_    Imaginary time.
     * @param itcfHF_ HF imaginary-time correlation function at this wave-vector
     *                and imaginary time tau.
     * @param lfc_    Span over the local field correction array (only lfc[0]
     *                is used).
     * @param itg_    Shared pointer to a 1D integrator.
     */
    ItcfGround(const double &x_,
               const std::shared_ptr<const Input> in_,
               const double &tau_,
               const double &itcfHF_,
               std::span<const double> lfc_,
               std::shared_ptr<Integrator1D> itg_);

    /** @brief Compute and return the ground-state RPA ITCF value. */
    double get();

  private:

    /** @brief Hartree-Fock contribution. */
    const double itcfHF;
    /** @brief Local field correction values. */
    std::span<const double> lfc;
    /** @brief 1D numerical integrator. */
    const std::shared_ptr<Integrator1D> itg;

    /**
     * @brief Integrand for the zero-temperature frequency integral.
     * @param Omega Real frequency value.
     * @return Value of the integrand at @p Omega.
     */
    double integrand(const double &Omega) const;
  };

} // namespace thermoUtil

#endif
