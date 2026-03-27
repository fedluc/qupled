#ifndef IET_HPP
#define IET_HPP

#include "schemes/input.hpp"
#include "schemes/stls.hpp"
#include "util/logger.hpp"
#include "util/numerics.hpp"
#include "util/vector2D.hpp"

/**
 * @brief Computes the bridge function for IET-based dielectric schemes.
 *
 * The IET (Integral Equation Theory) extension adds a bridge function to the
 * STLS local field correction. Supported bridge-function theories are HNC
 * (hypernetted-chain), IOI (Ichimaru–Ogata–Ichimaru), and LCT (Lucco
 * Castello–Tolias). The quantum-to-classical mapping for the effective coupling
 * parameter is selected via the input mapping string.
 */
class Iet : public Logger {

public:

  /**
   * @brief Construct the IET object.
   * @param in_   Shared pointer to the IET input parameters.
   * @param wvg_  Wave-vector grid inherited from the parent scheme.
   */
  explicit Iet(const std::shared_ptr<const IetInput> &in_,
               const std::vector<double> &wvg_)
      : Logger(true),
        inPtr(in_),
        wvg(wvg_),
        bf(wvg_.size()) {}

  /** @brief Compute the bridge function over the wave-vector grid. */
  void init();

  /** @brief Return the bridge function values over the wave-vector grid. */
  const std::vector<double> &getBf() const { return bf; }

  /**
   * @brief Attempt to load the initial LFC guess from the input parameters.
   * @param lfc Output LFC array to populate if a guess is available.
   * @return True if a valid guess was found, false otherwise.
   */
  bool initialGuessFromInput(Vector2D &lfc);

private:

  /** @brief Shared pointer to the IET input parameters. */
  const std::shared_ptr<const IetInput> inPtr;

  /** @brief Access the input as an @p IetInput reference. */
  const IetInput &in() const { return *inPtr; }

  /** @brief Cast the input to an @p IterationInput reference for RPA
   * parameters. */
  const IterationInput &inRpa() const {
    return *StlsUtil::dynamic_pointer_cast<IetInput, IterationInput>(inPtr);
  }

  /** @brief Wave-vector grid (shared with the parent scheme). */
  const std::vector<double> wvg;

  /** @brief Bridge function values over the wave-vector grid. */
  std::vector<double> bf;

  /** @brief Compute the bridge function for each wave-vector in @p wvg. */
  void computeBf();
};

/** @brief Internal helpers for the IET bridge function calculation. */
namespace IetUtil {

  /**
   * @brief Computes the bridge function at a single wave-vector.
   *
   * Supports three bridge-function theories (HNC, IOI, LCT) and three
   * quantum-to-classical mappings (coupling-parameter, energy, or
   * pressure-based). The effective classical coupling parameter is computed
   * first, then the corresponding bridge function is evaluated.
   */
  class BridgeFunction {

  public:

    /**
     * @brief Construct for a single wave-vector evaluation.
     * @param theory_  Name of the bridge-function theory ("HNC", "IOI", "LCT").
     * @param mapping_ Name of the quantum-classical mapping ("coupling", ...).
     * @param rs_      Quantum coupling parameter (Wigner–Seitz radius).
     * @param Theta_   Degeneracy parameter (reduced temperature).
     * @param x_       Dimensionless wave-vector.
     * @param itg_     Shared pointer to a 1D integrator (used by LCT).
     */
    BridgeFunction(const std::string &theory_,
                   const std::string &mapping_,
                   const double &rs_,
                   const double &Theta_,
                   const double &x_,
                   std::shared_ptr<Integrator1D> itg_)
        : theory(theory_),
          mapping(mapping_),
          rs(rs_),
          Theta(Theta_),
          x(x_),
          itg(itg_) {}

    /**
     * @brief Compute and return the bridge function.
     * @return Bridge-function value at @p x_.
     */
    double get() const;

  private:

    /** @brief Bridge-function theory identifier. */
    const std::string theory;
    /** @brief Quantum-classical mapping identifier. */
    const std::string mapping;
    /** @brief Quantum coupling parameter. */
    const double rs;
    /** @brief Degeneracy parameter. */
    const double Theta;
    /** @brief Wave-vector. */
    const double x;
    /** @brief 1D numerical integrator (used by the LCT theory). */
    const std::shared_ptr<Integrator1D> itg;
    /** @brief Unit-conversion constant \f$\lambda = (4/9\pi)^{1/3}\f$. */
    const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);

    /**
     * @brief Hypernetted-chain (HNC) bridge function.
     * @return HNC bridge-function value (zero by definition).
     */
    double hnc() const;

    /**
     * @brief Ichimaru–Ogata–Ichimaru (IOI) bridge function.
     * @return IOI bridge-function value at @p x_.
     */
    double ioi() const;

    /**
     * @brief Lucco Castello–Tolias (LCT) bridge function.
     * @return LCT bridge-function value at @p x_.
     */
    double lct() const;

    /**
     * @brief LCT integrand as a function of real-space distance @p r and
     *        classical coupling parameter @p Gamma.
     * @param r     Real-space distance.
     * @param Gamma Classical coupling parameter.
     * @return Integrand value.
     */
    double lctIntegrand(const double &r, const double &Gamma) const;

    /**
     * @brief Compute the effective classical coupling parameter for the bridge
     * function.
     * @return Classical coupling parameter \f$\Gamma\f$.
     */
    double couplingParameter() const;
  };

} // namespace IetUtil

#endif
