#ifndef ESA_HPP
#define ESA_HPP

#include "schemes/rpa.hpp"
#include <cmath>

// Forward declarations
class Dual22;

/**
 * @brief Solver for the ESA (Enhanced Screening Approximation) dielectric
 * scheme.
 *
 * Extends the RPA solver by replacing the zero local field correction with a
 * neural-network parametrization of the static local field correction (SLFC).
 * The parametrization satisfies the compressibility sum rule and uses
 * machine-learned coefficients fitted to quantum Monte Carlo data.
 */
class ESA : public Rpa {

public:

  /**
   * @brief Construct the ESA solver.
   * @param in_ Shared pointer to the input parameters.
   */
  explicit ESA(const std::shared_ptr<const Input> &in_)
      : Rpa(in_) {}

private:

  /** @brief Compute the neural-network-parametrized local field correction. */
  void computeLfc() override;
};

/** @brief Internal helpers for the ESA static local field correction. */
namespace ESAUtil {

  /**
   * @brief Computes the ESA static local field correction for a given state
   * point.
   *
   * The SLFC is constructed from a neural network parametrization that blends
   * the compressibility-sum-rule (CSR) result at long wavelengths with a
   * machine-learned fit at shorter wavelengths via a smooth activation
   * function.
   */
  class Slfc {

  public:

    /**
     * @brief Construct for a given coupling and degeneracy parameter.
     * @param rs_    Coupling parameter (Wigner–Seitz radius).
     * @param theta_ Degeneracy parameter (reduced temperature).
     */
    explicit Slfc(const double &rs_, const double &theta_)
        : rs(rs_),
          theta(theta_) {}

    /**
     * @brief Evaluate the SLFC at wave-vector @p x.
     * @param x Dimensionless wave-vector.
     * @return SLFC value at @p x.
     */
    double get(const double &x);

  public:

    /** @brief Coupling parameter. */
    const double rs;
    /** @brief Degeneracy parameter. */
    const double theta;

    /**
     * @brief Cache of pre-computed SLFC expansion coefficients.
     *
     * The coefficients are computed once and reused for all wave-vector
     * evaluations at the same state point.
     */
    struct Coefficients {
      /** @brief True if the coefficients are up to date and valid. */
      bool valid = false;
      /** @brief Long-wavelength limit coefficient. */
      double lwl;
      /** @brief Activation-function width parameter. */
      double afEta;
      /** @brief Activation-function inflection point. */
      double afxm;
      /** @brief Neural-network coefficient a. */
      double nna;
      /** @brief Neural-network coefficient b. */
      double nnb;
      /** @brief Neural-network coefficient c. */
      double nnc;
      /** @brief Neural-network coefficient d. */
      double nnd;
      /** @brief Compressibility-sum-rule coefficient. */
      double csr;
    };

    /** @brief Cached coefficient set for the current state point. */
    Coefficients coeff;

    /**
     * @brief Neural-network parametrization of the SLFC.
     * @param x Dimensionless wave-vector.
     * @return Neural-network contribution to the SLFC at @p x.
     */
    double nn(const double &x) const;

    /**
     * @brief Compressibility-sum-rule contribution to the SLFC.
     * @param x Dimensionless wave-vector.
     * @return CSR contribution at @p x.
     */
    double csr(const double &x) const;

    /** @brief Compute and cache all SLFC coefficients. */
    void computeCoefficients();

    /** @brief Compute and cache the neural-network coefficients. */
    void computeNNCoefficients();

    /** @brief Compute and cache the compressibility-sum-rule coefficient. */
    void computeCSRCoefficients();

    /**
     * @brief On-top value of the radial distribution function.
     * @return g(0) at the current state point.
     */
    double onTop() const;

    /**
     * @brief Activation function that blends the long- and short-wavelength
     * limits.
     * @param x Dimensionless wave-vector.
     * @return Value of the activation function at @p x.
     */
    double activationFunction(const double &x) const;

    /**
     * @brief Evaluate the parametrized free energy (value only).
     * @return Free energy dual number at the current (rs, theta).
     */
    Dual22 freeEnergy() const;

    /**
     * @brief Evaluate the parametrized free energy at arbitrary (rs, theta).
     * @param rs    Dual-number coupling parameter.
     * @param theta Dual-number degeneracy parameter.
     * @return Free energy dual number carrying value and derivatives.
     */
    Dual22 freeEnergy(const Dual22 &rs, const Dual22 &theta) const;
  };

} // namespace ESAUtil
#endif
