#ifndef STLSIET_HPP
#define STLSIET_HPP

#include "schemes/iet.hpp"
#include "schemes/stls.hpp"

/**
 * @brief Solver for the STLS-IET dielectric scheme.
 *
 * Combines the STLS iterative solver with an IET bridge function. The static
 * local field correction includes both the STLS exchange-correlation term and
 * a bridge-function contribution computed by the @p Iet helper. Convergence
 * is driven by the same iterative mixing as for plain STLS.
 */
class StlsIet : public Stls {

public:

  /**
   * @brief Construct the STLS-IET solver.
   * @param in_ Shared pointer to the STLS-IET input parameters.
   */
  explicit StlsIet(const std::shared_ptr<const StlsIetInput> &in_);

  /**
   * @brief Return the bridge function values over the wave-vector grid.
   * @return Constant reference to the bridge function array.
   */
  const std::vector<double> &getBf() const { return iet.getBf(); }

private:

  /** @brief IET helper that computes the bridge function. */
  Iet iet;
  /** @brief 2D integrator for the SLFC double integral. */
  const std::shared_ptr<Integrator2D> itg2D;
  /** @brief Quadrature grid for the 2D integral. */
  std::vector<double> itgGrid;

  /** @brief Access the input as a @p StlsIetInput reference. */
  const StlsIetInput &in() const {
    return *StlsUtil::dynamic_pointer_cast<Input, StlsIetInput>(inPtr);
  }

  /** @brief Initialize wave-vector grid, chemical potential, and bridge
   * function. */
  void init() override;

  /** @brief Compute the STLS-IET static local field correction. */
  void computeLfc() override;

  /**
   * @brief Attempt to load the initial LFC guess from the input parameters.
   * @return True if a valid guess was found, false otherwise.
   */
  bool initialGuessFromInput() override;
};

/** @brief Internal helpers for the STLS-IET static local field correction. */
namespace StlsIetUtil {

  /**
   * @brief Computes the STLS-IET static local field correction.
   *
   * Evaluates the SLFC integral that combines the standard STLS exchange-
   * correlation term with a bridge-function contribution via a 2D quadrature.
   */
  class Slfc : public StlsUtil::SlfcBase, dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a STLS-IET SLFC calculation.
     * @param x_        Wave-vector value.
     * @param yMin_     Lower integration limit.
     * @param yMax_     Upper integration limit.
     * @param ssfi_     Shared pointer to an SSF interpolator.
     * @param lfci_     Shared pointer to an LFC interpolator.
     * @param bfi_      Shared pointer to a bridge-function interpolator.
     * @param itgGrid_  Grid for 2D integration.
     * @param itg_      Shared pointer to a 2D integrator.
     * @param in_       Shared pointer to the input parameters.
     */
    Slfc(const double &x_,
         const double &yMin_,
         const double &yMax_,
         std::shared_ptr<Interpolator1D> ssfi_,
         std::shared_ptr<Interpolator1D> lfci_,
         std::shared_ptr<Interpolator1D> bfi_,
         const std::vector<double> &itgGrid_,
         std::shared_ptr<Integrator2D> itg_,
         const std::shared_ptr<const Input> in_)
        : SlfcBase(x_, yMin_, yMax_, ssfi_),
          itg(itg_),
          itgGrid(itgGrid_),
          lfci(lfci_),
          bfi(bfi_),
          res(x_),
          in(in_) {}

    /**
     * @brief Compute and return the SLFC at the current wave-vector.
     * @return SLFC value including the bridge-function contribution.
     */
    double get();

  private:

    /** @brief 2D numerical integrator. */
    const std::shared_ptr<Integrator2D> itg;
    /** @brief Grid for 2D integration. */
    const std::vector<double> itgGrid;
    /**
     * @brief 3D outer integrand over auxiliary momentum @p y.
     * @param y Outer integration variable.
     */
    double integrand1(const double &y) const;
    /**
     * @brief 3D inner integrand over auxiliary momentum @p w.
     * @param w Inner integration variable.
     */
    double integrand2(const double &w) const;
    /**
     * @brief 2D outer integrand over auxiliary momentum @p y.
     * @param y Outer integration variable.
     */
    double integrand1_2D(const double &y) const;
    /**
     * @brief 2D inner integrand over auxiliary momentum @p w.
     * @param w Inner integration variable.
     */
    double integrand2_2D(const double &w) const;
    /** @brief Interpolator for the static local field correction. */
    const std::shared_ptr<Interpolator1D> lfci;
    /** @brief Interpolator for the bridge function. */
    const std::shared_ptr<Interpolator1D> bfi;
    /** @brief Result of the SLFC integration. */
    double res;
    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;
    /**
     * @brief Evaluate the interpolated LFC at wave-vector @p x.
     * @param x Wave-vector value.
     * @return Interpolated LFC value.
     */
    double lfc(const double &x) const;
    /**
     * @brief Evaluate the interpolated bridge function at wave-vector @p x_.
     * @param x_ Wave-vector value.
     * @return Interpolated bridge function value.
     */
    double bf(const double &x_) const;
    void compute2D() override;
    void compute3D() override;
  };

} // namespace StlsIetUtil

#endif
