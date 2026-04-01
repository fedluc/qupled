#ifndef RPA_HPP
#define RPA_HPP

#include "schemes/hf.hpp"
#include "schemes/input.hpp"
#include "util/logger.hpp"
#include "util/numerics.hpp"
#include "util/vector2D.hpp"
#include <vector>

/**
 * @brief Solver for the Random Phase Approximation (RPA) dielectric scheme.
 *
 * Extends the Hartree-Fock (HF) base class to compute the static structure
 * factor (SSF) and the imaginary-time correlation function (ITCF) within the
 * RPA. The SSF is obtained as the special case ITCF(tau=0). Finite-temperature
 * and zero-temperature (ground state) regimes are both supported; the ITCF is
 * only implemented for 3D finite-temperature systems.
 */
class Rpa : public HF {

public:

  /**
   * @brief Construct with explicit verbosity flag.
   * @param in_       Shared pointer to the input parameters.
   * @param verbose_  If false, solver output is suppressed.
   */
  Rpa(const std::shared_ptr<const Input> &in_, const bool verbose_);

  /**
   * @brief Construct with verbosity enabled.
   * @param in_ Shared pointer to the input parameters.
   */
  explicit Rpa(const std::shared_ptr<const Input> &in_)
      : Rpa(in_, true) {}

protected:

  /** @brief Hartree-Fock static structure factor. */
  std::vector<double> ssfHF;

  /** @brief Initialize wave-vector grid and other basic properties. */
  void init() override;

  /** @brief Compute the static structure factor at finite temperature. */
  void computeSsfFinite() override;

  /** @brief Compute the static structure factor at zero temperature. */
  void computeSsfGround() override;

private:

  /** @brief Compute the Hartree-Fock static structure factor. */
  void computeSsfHF();

  /** @brief Compute the local field correction (zero for RPA). */
  void computeLfc() override;
};

/** @brief Internal helpers for the RPA static structure factor computation. */
namespace RpaUtil {

  /**
   * @brief Base class holding shared state for static structure factor helpers.
   */
  class SsfBase {

  protected:

    /**
     * @brief Construct the base with the quantities needed for SSF evaluation.
     * @param x_     Wave-vector value.
     * @param ssfHF_ Hartree-Fock static structure factor at this wave-vector.
     * @param lfc_   Span over the local field correction array.
     * @param in_    Shared pointer to the input parameters.
     */
    SsfBase(const double &x_,
            const double &ssfHF_,
            std::span<const double> lfc_,
            const std::shared_ptr<const Input> in_)
        : x(x_),
          ssfHF(ssfHF_),
          lfc(lfc_),
          in(in_) {}

    /** @brief Wave-vector value. */
    const double x;

    /** @brief Hartree-Fock contribution to the static structure factor. */
    const double ssfHF;

    /** @brief Local field correction values. */
    std::span<const double> lfc;

    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;

    /** @brief Normalized interaction potential at the current wave-vector. */
    double ip() const;
  };

  /**
   * @brief Computes the finite-temperature RPA static structure factor.
   *
   * Integrates over Matsubara frequencies using the ideal density response.
   */
  class Ssf : public SsfBase, dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a finite-temperature calculation.
     * @param x_     Wave-vector value.
     * @param ssfHF_ Hartree-Fock static structure factor at this wave-vector.
     * @param lfc_   Span over the local field correction array.
     * @param in_    Shared pointer to the input parameters.
     * @param idr_   Span over the ideal density response array.
     */
    Ssf(const double &x_,
        const double &ssfHF_,
        std::span<const double> lfc_,
        const std::shared_ptr<const Input> in_,
        std::span<const double> idr_)
        : SsfBase(x_, ssfHF_, lfc_, in_),
          idr(idr_),
          res(numUtil::NaN) {}

    /** @brief Compute and return the static structure factor. */
    double get();

  protected:

    /** @brief Ideal density response values over Matsubara frequencies. */
    const std::span<const double> idr;

  private:

    /** @brief Stores the result of the frequency summation. */
    double res;

    void compute2D() override;
    void compute3D() override;
  };

  /**
   * @brief Computes the finite-temperature RPA imaginary-time correlation
   * function (ITCF).
   *
   * Evaluates F(x, tau) = F_HF(x, tau) minus the Matsubara correction sum
   * weighted by cos(2*pi*l*tau). The SSF is recovered as the special case
   * tau = 0 by delegating to the Ssf class.
   */
  class Itcf : public SsfBase, dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a finite-temperature ITCF calculation.
     * @param x_      Wave-vector value.
     * @param itcfHF_ HF imaginary-time correlation function at this
     * wave-vector.
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
         const double &tau_)
        : SsfBase(x_, itcfHF_, lfc_, in_),
          idr(idr_),
          tau(tau_),
          res(numUtil::NaN) {}

    /** @brief Compute and return the RPA ITCF value. */
    double get();

  private:

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

  /**
   * @brief Computes the zero-temperature RPA static structure factor.
   *
   * Uses a 1D numerical integrator over real frequencies.
   */
  class SsfGround : public SsfBase {

  public:

    /**
     * @brief Construct for a zero-temperature calculation.
     * @param x_     Wave-vector value.
     * @param ssfHF_ Hartree-Fock static structure factor at this wave-vector.
     * @param lfc_   Span over the local field correction array.
     * @param itg_   Shared pointer to a 1D integrator instance.
     * @param in_    Shared pointer to the input parameters.
     */
    SsfGround(const double &x_,
              const double &ssfHF_,
              std::span<const double> lfc_,
              std::shared_ptr<Integrator1D> itg_,
              const std::shared_ptr<const Input> in_)
        : SsfBase(x_, ssfHF_, lfc_, in_),
          itg(itg_) {}

    /** @brief Compute and return the static structure factor. */
    double get();

  protected:

    /** @brief 1D numerical integrator. */
    const std::shared_ptr<Integrator1D> itg;

    /**
     * @brief Integrand for the zero-temperature frequency integral.
     * @param Omega Real frequency value.
     * @return Value of the integrand at @p Omega.
     */
    double integrand(const double &Omega) const;
  };

} // namespace RpaUtil

#endif
