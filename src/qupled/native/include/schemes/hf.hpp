#ifndef HF_HPP
#define HF_HPP

#include "schemes/input.hpp"
#include "util/dimensions_util.hpp"
#include "util/logger.hpp"
#include "util/numerics.hpp"
#include <vector>

/**
 * @brief Solver for the Hartree-Fock (HF) dielectric scheme.
 *
 * Base class for all dielectric scheme solvers. Computes the ideal density
 * response (IDR), the static structure factor (SSF), and the local field
 * correction (LFC) on a wave-vector grid. Derived classes override the
 * virtual compute methods to implement higher-level approximations.
 */
class HF : public Logger {

public:

  /**
   * @brief Construct with explicit verbosity flag.
   * @param in_       Shared pointer to the input parameters.
   * @param verbose_  If false, solver output is suppressed.
   */
  HF(const std::shared_ptr<const Input> &in_, const bool verbose_);

  /**
   * @brief Construct with verbosity enabled.
   * @param in_ Shared pointer to the input parameters.
   */
  explicit HF(const std::shared_ptr<const Input> &in_)
      : HF(in_, true) {}

  /** @brief Virtual destructor. */
  virtual ~HF() = default;

  /**
   * @brief Run the full solver pipeline.
   * @return 0 on success.
   */
  int compute();

  /** @brief Return the ideal density response grid (wave-vectors × Matsubara
   * frequencies). */
  const Vector2D &getIdr() const { return idr; }

  /** @brief Return the local field correction grid (wave-vectors × Matsubara
   * frequencies). */
  const Vector2D &getLfc() const { return lfc; }

  /** @brief Return the static structure factor over the wave-vector grid. */
  const std::vector<double> &getSsf() const { return ssf; }

  /** @brief Return the wave-vector grid. */
  const std::vector<double> &getWvg() const { return wvg; }

  /**
   * @brief Return the chemical potential.
   * @return Chemical potential (in units of the thermal energy).
   */
  double getChemicalPotential() const { return mu; }

  /**
   * @brief Compute and return the static density response.
   * @return Vector of static density response values over the wave-vector grid.
   */
  std::vector<double> getSdr() const;

  /**
   * @brief Compute and return the interaction energy.
   * @return Dimensionless interaction energy per particle.
   */
  double getUInt() const;

protected:

  /** @brief Shared pointer to the input parameters. */
  const std::shared_ptr<const Input> inPtr;

  /** @brief 1D numerical integrator. */
  const std::shared_ptr<Integrator1D> itg;

  /** @brief Wave-vector grid. */
  std::vector<double> wvg;

  /** @brief Ideal density response (rows = wave-vectors, columns = Matsubara
   * frequencies). */
  Vector2D idr;

  /** @brief Static local field correction (rows = wave-vectors, columns =
   * Matsubara frequencies). */
  Vector2D lfc;

  /** @brief Static structure factor over the wave-vector grid. */
  std::vector<double> ssf;

  /** @brief Chemical potential (in units of the thermal energy). */
  double mu;

  /** @brief Dereference the input shared pointer. */
  const Input &in() const { return *inPtr; }

  /** @brief Initialize the wave-vector grid and other basic derived quantities.
   */
  virtual void init();

  /** @brief Compute all structural properties (IDR, SSF, LFC). */
  virtual void computeStructuralProperties();

  /** @brief Dispatch the SSF calculation to the appropriate temperature branch.
   */
  virtual void computeSsf();

  /** @brief Compute the static structure factor at finite temperature. */
  virtual void computeSsfFinite();

  /** @brief Compute the static structure factor at zero temperature (ground
   * state). */
  virtual void computeSsfGround();

protected:

  /** @brief Compute the local field correction (zero for bare HF). */
  virtual void computeLfc();

private:

  /** @brief Construct the uniform wave-vector grid from the input parameters.
   */
  void buildWaveVectorGrid();

  /** @brief Compute the chemical potential from the normalization condition. */
  void computeChemicalPotential();

  /** @brief Dispatch the IDR calculation to the appropriate temperature branch.
   */
  void computeIdr();

  /** @brief Compute the ideal density response at finite temperature. */
  void computeIdrFinite();

  /** @brief Compute the ideal density response at zero temperature. */
  void computeIdrGround();
};

/** @brief Internal helpers for the Hartree-Fock ideal density response and SSF.
 */
namespace HFUtil {

  /**
   * @brief Computes the ideal density response at finite temperature.
   *
   * Evaluates the IDR for a given wave-vector @p x_ over all Matsubara
   * frequencies by numerical integration over the auxiliary momentum.
   */
  class Idr : public dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a finite-temperature IDR calculation.
     * @param in_    Shared pointer to the input parameters.
     * @param x_     Wave-vector value.
     * @param mu_    Chemical potential.
     * @param yMin_  Lower integration limit.
     * @param yMax_  Upper integration limit.
     * @param itg_   Shared pointer to a 1D integrator.
     */
    Idr(const std::shared_ptr<const Input> in_,
        const double &x_,
        const double &mu_,
        const double &yMin_,
        const double &yMax_,
        std::shared_ptr<Integrator1D> itg_)
        : in(in_),
          x(x_),
          mu(mu_),
          yMin(yMin_),
          yMax(yMax_),
          itg(itg_),
          res(in_->getNMatsubara()) {}

    /**
     * @brief Compute and return the IDR for all Matsubara frequencies.
     * @return Vector of IDR values, one per Matsubara frequency.
     */
    std::vector<double> get();

  private:

    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;
    /** @brief Wave-vector. */
    const double x;
    /** @brief Chemical potential. */
    const double mu;
    /** @brief Lower integration limit. */
    const double yMin;
    /** @brief Upper integration limit. */
    const double yMax;
    /** @brief 1D numerical integrator. */
    const std::shared_ptr<Integrator1D> itg;
    /** @brief Storage for the per-frequency integration results. */
    std::vector<double> res;
    void compute3D() override;
    void compute2D() override;
    /**
     * @brief 3D integrand for Matsubara frequency @p l.
     * @param y Auxiliary momentum variable.
     * @param l Matsubara frequency index.
     */
    double integrand(const double &y, const int &l) const;
    /**
     * @brief 3D integrand for the zeroth Matsubara frequency.
     * @param y Auxiliary momentum variable.
     */
    double integrand(const double &y) const;
    /**
     * @brief 2D integrand for Matsubara frequency @p l.
     * @param y Auxiliary momentum variable.
     * @param l Matsubara frequency index.
     */
    double integrand2D(const double &y, const int &l) const;
    /**
     * @brief 2D integrand for the zeroth Matsubara frequency.
     * @param y Auxiliary momentum variable.
     */
    double integrand2D(const double &y) const;
  };

  /**
   * @brief Computes the ideal density response at zero temperature.
   *
   * Uses the analytic zero-temperature Lindhard function evaluated at
   * real frequency @p Omega_ and wave-vector @p x_.
   */
  class IdrGround : public dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a ground-state IDR calculation.
     * @param x_     Wave-vector value.
     * @param Omega_ Real frequency value.
     */
    IdrGround(const std::shared_ptr<const Input> in_,
              const double &x_,
              const double &Omega_)
        : in(in_),
          x(x_),
          Omega(Omega_),
          itg(std::make_shared<Integrator1D>(Integrator1D::Type::DEFAULT,
                                             1e-6)) {}

    /**
     * @brief Compute and return the ground-state IDR.
     * @return Analytic IDR value at (@p x_, @p Omega_).
     */
    double get();

  private:

    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;
    /** @brief Wave-vector. */
    const double x;
    /** @brief Real frequency. */
    const double Omega;
    /** @brief 1D numerical integrator. */
    const std::shared_ptr<Integrator1D> itg;
    /** @brief Result of the ground-state IDR computation. */
    double res;
    void compute3D() override;
    void compute2D() override;
    /**
     * @brief 2D integrand.
     * @param y Auxiliary momentum variable.
     */
    double integrand2D(const double &y) const;
  };

  /**
   * @brief Computes the Hartree-Fock static structure factor at finite
   * temperature.
   *
   * Evaluates F_HF(x, 0) via numerical integration over the auxiliary
   * momentum. For 3D systems, uses a single 1D integral. For 2D systems,
   * uses a 2D integral (over y and angle p) plus the IDR contribution.
   */
  class Ssf : public dimensionsUtil::DimensionsHandler {

  public:

    /**
     * @brief Construct for a finite-temperature SSF calculation.
     * @param in_        Shared pointer to the input parameters.
     * @param x_         Wave-vector value.
     * @param mu_        Chemical potential.
     * @param yMin_      Lower integration limit.
     * @param yMax_      Upper integration limit.
     * @param itg_       Shared pointer to a 1D integrator.
     * @param itgGrid_   Grid for 2D integration.
     * @param itg2_      Shared pointer to a 2D integrator.
     * @param idr0_      Ideal density response at l=0 for the current
     * wave-vector.
     */
    Ssf(const std::shared_ptr<const Input> in_,
        const double &x_,
        const double &mu_,
        const double &yMin_,
        const double &yMax_,
        std::shared_ptr<Integrator1D> itg_,
        const std::vector<double> &itgGrid_,
        std::shared_ptr<Integrator2D> itg2_,
        const double &idr0_)
        : in(in_),
          x(x_),
          mu(mu_),
          yMin(yMin_),
          yMax(yMax_),
          itg(itg_),
          itgGrid(itgGrid_),
          itg2(itg2_),
          idr0(idr0_),
          res(x_) {}

    /**
     * @brief Compute and return the HF static structure factor.
     * @return SSF value at the current wave-vector.
     */
    double get();

  private:

    /** @brief Input parameters. */
    const std::shared_ptr<const Input> in;
    /** @brief Wave-vector. */
    const double x;
    /** @brief Chemical potential. */
    const double mu;
    /** @brief Lower integration limit. */
    const double yMin;
    /** @brief Upper integration limit. */
    const double yMax;
    /** @brief 1D numerical integrator. */
    const std::shared_ptr<Integrator1D> itg;
    /** @brief Grid for 2D integration. */
    const std::vector<double> &itgGrid;
    /** @brief 2D numerical integrator. */
    const std::shared_ptr<Integrator2D> itg2;
    /** @brief Ideal density response at l=0 for the current wave-vector. */
    const double idr0;
    /** @brief Result of the SSF computation. */
    double res;
    void compute3D() override;
    void compute2D() override;
    /**
     * @brief 3D integrand over auxiliary momentum @p y.
     * @param y Auxiliary momentum variable.
     */
    double integrand(const double &y) const;
    /**
     * @brief Outer 2D integrand over auxiliary momentum @p y.
     * @param y Outer integration variable.
     */
    double integrand2DOut(const double &y) const;
    /**
     * @brief Inner 2D integrand (angle) over @p p.
     * @param p Inner integration variable (angle).
     */
    double integrand2DIn(const double &p) const;
  };

  /**
   * @brief Computes the Hartree-Fock static structure factor at zero
   * temperature.
   *
   * Uses the analytic ground-state expression for the HF SSF.
   */
  class SsfGround {

  public:

    /**
     * @brief Construct for a zero-temperature SSF calculation.
     * @param x_ Wave-vector value.
     */
    explicit SsfGround(const double &x_)
        : x(x_) {}

    /**
     * @brief Compute and return the ground-state HF static structure factor.
     * @return Analytic SSF value at @p x_.
     */
    double get() const;

  private:

    /** @brief Wave-vector. */
    const double x;
  };

} // namespace HFUtil

#endif
