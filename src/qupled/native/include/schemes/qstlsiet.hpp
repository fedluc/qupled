#ifndef QSTLSIET_HPP
#define QSTLSIET_HPP

#include "schemes/iet.hpp"
#include "schemes/input.hpp"
#include "schemes/qstls.hpp"

/**
 * @brief Solver for the qSTLS-IET (quantum STLS with Integral Equation Theory)
 * scheme.
 *
 * Extends the qSTLS solver by incorporating an IET bridge-function contribution
 * into the local field correction. Both the quantum (ADR-derived) and the
 * classical (bridge-function) parts of the LFC are iterated to
 * self-consistency.
 */
class QstlsIet : public Qstls {

public:

  /**
   * @brief Construct the qSTLS-IET solver.
   * @param in_ Shared pointer to the qSTLS-IET input parameters.
   */
  explicit QstlsIet(const std::shared_ptr<const QstlsIetInput> &in_);

  /**
   * @brief Return the bridge function values over the wave-vector grid.
   * @return Constant reference to the bridge function array.
   */
  const std::vector<double> &getBf() const { return iet.getBf(); }

private:

  /** @brief IET helper that computes the bridge function. */
  Iet iet;
  /** @brief IET contribution to the local field correction. */
  Vector2D lfcIet;
  /** @brief Quadrature grid for the 2D ADR integral. */
  std::vector<double> itgGrid;

  /** @brief Access the input as a @p QstlsIetInput reference. */
  const QstlsIetInput &in() const {
    return *StlsUtil::dynamic_pointer_cast<Input, QstlsIetInput>(inPtr);
  }

  /** @brief Initialize wave-vector grid, chemical potential, bridge function,
   * and fixed ADR. */
  void init() override;

  /** @brief Compute the qSTLS-IET local field correction. */
  void computeLfc() override;

  /** @brief Compute the fixed (frequency-independent) ADR component. */
  void computeAdrFixed();

  /**
   * @brief Attempt to load the initial LFC guess from the input parameters.
   * @return True if a valid guess was found, false otherwise.
   */
  bool initialGuessFromInput() override;
};

/** @brief Internal helpers for the qSTLS-IET auxiliary density response. */
namespace QstlsIetUtil {

  /**
   * @brief Computes the IET contribution to the finite-temperature ADR.
   *
   * Extends the standard qSTLS ADR by adding bridge-function and dynamic LFC
   * terms. The integration combines the fixed ADR component with the current
   * SSF, LFC, and bridge function.
   */
  class AdrIet : public QstlsUtil::AdrBase {

  public:

    /**
     * @brief Construct for a finite-temperature IET-ADR calculation.
     * @param Theta_    Degeneracy parameter.
     * @param qMin_     Lower integration limit.
     * @param qMax_     Upper integration limit.
     * @param x_        Wave-vector value.
     * @param ssfi_     Shared pointer to an SSF interpolator.
     * @param lfci_     Per-frequency LFC interpolators (one per Matsubara
     * frequency).
     * @param bfi_      Shared pointer to a bridge-function interpolator.
     * @param itgGrid_  Grid for 2D integration.
     * @param itg_      Shared pointer to a 2D integrator.
     */
    AdrIet(const double &Theta_,
           const double &qMin_,
           const double &qMax_,
           const double &x_,
           std::shared_ptr<Interpolator1D> ssfi_,
           std::vector<std::shared_ptr<Interpolator1D>> lfci_,
           std::shared_ptr<Interpolator1D> bfi_,
           const std::vector<double> &itgGrid_,
           std::shared_ptr<Integrator2D> itg_)
        : QstlsUtil::AdrBase(Theta_, qMin_, qMax_, x_, ssfi_),
          itg(itg_),
          itgGrid(itgGrid_),
          lfci(lfci_),
          bfi(bfi_) {}

    /**
     * @brief Compute and store the IET-ADR for all Matsubara frequencies.
     * @param wvg   Wave-vector grid.
     * @param fixed Fixed ADR component (wave-vectors × wave-vectors ×
     * frequencies).
     * @param res   Output array (wave-vectors × frequencies) to accumulate
     * into.
     */
    void
    get(const std::vector<double> &wvg, const Vector3D &fixed, Vector2D &res);

  private:

    /** @brief Reference to lower integration limit (alias for @p yMin). */
    const double &qMin = yMin;
    /** @brief Reference to upper integration limit (alias for @p yMax). */
    const double &qMax = yMax;

    /**
     * @brief Outer integrand over auxiliary momentum @p q and Matsubara index
     * @p l.
     * @param q Auxiliary momentum.
     * @param l Matsubara frequency index.
     */
    double integrand1(const double &q, const int &l) const;

    /**
     * @brief Inner integrand over auxiliary momentum @p y.
     * @param y Auxiliary momentum.
     */
    double integrand2(const double &y) const;

    /** @brief 2D numerical integrator. */
    const std::shared_ptr<Integrator2D> itg;
    /** @brief Grid for 2D integration. */
    const std::vector<double> &itgGrid;
    /** @brief Per-frequency LFC interpolators. */
    const std::vector<std::shared_ptr<Interpolator1D>> lfci;
    /** @brief Bridge function interpolator. */
    const std::shared_ptr<Interpolator1D> bfi;
    /** @brief 2D interpolator for the fixed ADR component. */
    Interpolator2D fixi;

    /**
     * @brief Evaluate the dynamic LFC at wave-vector @p y and frequency index
     * @p l.
     * @param y Auxiliary momentum.
     * @param l Matsubara frequency index.
     */
    double lfc(const double &y, const int &l) const;

    /**
     * @brief Evaluate the bridge-function contribution at wave-vector @p y.
     * @param y Wave-vector value.
     */
    double bf(const double &y) const;

    /**
     * @brief Evaluate the fixed ADR component at (@p x, @p y).
     * @param x Primary wave-vector.
     * @param y Auxiliary wave-vector.
     */
    double fix(const double &x, const double &y) const;
  };

  /**
   * @brief Computes the IET fixed component of the ADR.
   *
   * This component depends only on Fermi–Dirac occupation numbers and the
   * Matsubara frequency index, not on the self-consistent SSF or LFC.
   */
  class AdrFixedIet : public QstlsUtil::AdrFixedBase {

  public:

    /**
     * @brief Construct for a finite-temperature fixed-IET-ADR calculation.
     * @param Theta_ Degeneracy parameter.
     * @param qMin_  Lower integration limit.
     * @param qMax_  Upper integration limit.
     * @param x_     Wave-vector value.
     * @param mu_    Chemical potential.
     * @param itg_   Shared pointer to a 1D integrator.
     */
    AdrFixedIet(const double &Theta_,
                const double &qMin_,
                const double &qMax_,
                const double &x_,
                const double &mu_,
                std::shared_ptr<Integrator1D> itg_)
        : QstlsUtil::AdrFixedBase(Theta_, qMin_, qMax_, x_, mu_),
          itg(itg_) {}

    /**
     * @brief Compute and store the fixed IET-ADR for a single Matsubara
     * frequency.
     * @param l    Matsubara frequency index.
     * @param wvg  Wave-vector grid.
     * @param res  Output 3D array to accumulate into.
     */
    void get(int l, const std::vector<double> &wvg, Vector3D &res) const;

  private:

    /** @brief Reference to lower integration limit (alias for @p qMin). */
    const double &tMin = qMin;
    /** @brief Reference to upper integration limit (alias for @p qMax). */
    const double &tMax = qMax;

    /**
     * @brief Integrand over variables @p t, @p y, @p q, and @p l.
     * @param t Intermediate momentum variable.
     * @param y Auxiliary momentum variable.
     * @param q Auxiliary momentum variable.
     * @param l Matsubara frequency value.
     */
    double integrand(const double &t,
                     const double &y,
                     const double &q,
                     const double &l) const;

    /** @brief 1D numerical integrator. */
    const std::shared_ptr<Integrator1D> itg;
  };

} // namespace QstlsIetUtil

#endif
