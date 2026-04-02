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
 * factor (SSF) within the RPA. Both finite-temperature and zero-temperature
 * (ground state) regimes are supported. Helper classes for the SSF
 * computation are provided in the RpaUtil namespace.
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
   * @brief Base class holding shared state for SSF helpers.
   */
  class SsfBase {

  protected:

    /**
     * @brief Construct the base with the quantities needed for SSF evaluation.
     * @param x_     Wave-vector value.
     * @param ssfHF_ Hartree-Fock SSF at this wave-vector.
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

    /** @brief Hartree-Fock SSF. */
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
    /**
     * @brief Compute the Matsubara frequency summation for SSF.
     * @return Sum over Matsubara frequencies (unweighted, for tau=0).
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
