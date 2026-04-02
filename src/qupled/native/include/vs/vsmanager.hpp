#ifndef VS_VSMANAGER_HPP
#define VS_VSMANAGER_HPP

#include "schemes/input.hpp"
#include "vs/grid_point.hpp"
#include "vs/vsworker.hpp"
#include <array>
#include <memory>
#include <vector>

/**
 * @brief Manages the 3×3 grid of workers for the variational-swarm scheme.
 *
 * @p VSManager owns nine @p VSWorker instances arranged on a 3×3 grid in the
 * (rs, Theta) plane. It coordinates the self-consistency iterations by
 * driving each worker through the init / guess / SSF / LFC cycle and then
 * applying finite-difference LFC derivatives across the grid to implement the
 * variational-swarm local-field-correction correction.
 *
 * Concrete subclasses supply input accessors and, via @p compute(), the
 * per-scheme initialization of the @p workers array.
 */
class VSManager {
public:

  /**
   * @brief Set the free variational parameter alpha.
   * @param alpha_ New value of alpha.
   */
  void setAlpha(double alpha_) { alpha = alpha_; }

  /**
   * @brief Return the current value of the free parameter alpha.
   * @return Alpha value.
   */
  double getAlpha() const { return alpha; }

  /**
   * @brief Return the convergence error from the last iteration.
   * @return RMS difference between successive SSF iterates (maximum over all
   * workers).
   */
  double getError() const { return computeError(); }

  /**
   * @brief Run one full self-consistency sweep over all grid points.
   *
   * Implemented by derived classes to initialize workers and drive the
   * iterative loop until convergence.
   *
   * @return 0 on success.
   */
  virtual int compute() = 0;

  /**
   * @brief Return the SSF at a given grid point.
   * @param p Grid point identifier.
   * @return Constant reference to the SSF vector.
   */
  const std::vector<double> &getSsf(GridPoint p) const;

  /**
   * @brief Return the LFC at a given grid point.
   * @param p Grid point identifier.
   * @return Constant reference to the LFC array.
   */
  const Vector2D &getLfc(GridPoint p) const;

  /**
   * @brief Return the wave-vector grid at a given grid point.
   * @param p Grid point identifier.
   * @return Constant reference to the wave-vector grid.
   */
  const std::vector<double> &getWvg(GridPoint p) const;

  /**
   * @brief Return the ideal density response at a given grid point.
   * @param p Grid point identifier.
   * @return Constant reference to the IDR array.
   */
  const Vector2D &getIdr(GridPoint p) const;

  /**
   * @brief Compute and return the static density response at a given grid
   * point.
   * @param p Grid point identifier.
   * @return Static density response vector.
   */
  std::vector<double> getSdr(GridPoint p) const;

  /**
   * @brief Return the coupling parameter at a given grid point.
   * @param p Grid point identifier.
   * @return rs value.
   */
  double getCoupling(GridPoint p) const;

  /**
   * @brief Return the degeneracy parameter at a given grid point.
   * @param p Grid point identifier.
   * @return Theta value.
   */
  double getDegeneracy(GridPoint p) const;

  /**
   * @brief Return the interaction energy at a given grid point.
   * @param p Grid point identifier.
   * @return Dimensionless interaction energy per particle.
   */
  double getUInt(GridPoint p) const;

  /**
   * @brief Return the chemical potential at a given grid point.
   * @param p Grid point identifier.
   * @return Chemical potential (in units of the thermal energy).
   */
  double getChemicalPotential(GridPoint p) const;

  /**
   * @brief Return the free-energy integrand contribution (Q-adder) at a grid
   * point.
   * @param p Grid point identifier.
   * @return Q-adder value.
   */
  double getQAdder(GridPoint p) const;

  /**
   * @brief Return the free-energy integrand value at a given grid point.
   * @param p Grid point identifier.
   * @return Free-energy integrand value.
   */
  double getFxcIntegrandValue(GridPoint p) const;

protected:

  /** @brief Construct with alpha initialized to infinity and @p initDone =
   * false. */
  explicit VSManager()
      : alpha(numUtil::Inf),
        initDone(false) {}

  /**
   * @brief Metadata for finite-difference derivative computation at one grid
   * point.
   */
  struct DerivativeData {
    /** @brief Derivative type for the finite-difference stencil. */
    enum class Type { CENTERED, FORWARD, BACKWARD };
    /** @brief Stencil type to use at this grid point. */
    Type type;
    /** @brief Flat index of the "up" neighbour in the stencil. */
    size_t upIdx;
    /** @brief Flat index of the "down" neighbour in the stencil. */
    size_t downIdx;
  };

  /** @brief Number of grid points (3 × 3 = 9). */
  static constexpr int N = 9;

  /** @brief The nine worker objects, one per grid point. */
  std::array<std::unique_ptr<VSWorker>, N> workers;

  /** @brief LFC finite-difference derivatives, one per grid point. */
  std::array<Vector2D, N> lfcDerivatives;

  /** @brief Finite-difference metadata along the rs axis. */
  std::array<DerivativeData, N> rsDerivData;

  /** @brief Finite-difference metadata along the Theta axis. */
  std::array<DerivativeData, N> thetaDerivData;

  /** @brief rs values for each of the nine grid points. */
  std::array<double, N> rsValues;

  /** @brief Theta values for each of the nine grid points. */
  std::array<double, N> thetaValues;

  /** @brief Current value of the free variational parameter alpha. */
  double alpha;

  /** @brief True after the first call to @p init(); prevents re-initialization.
   */
  bool initDone;

  /** @brief Populate @p rsDerivData and @p thetaDerivData from the grid layout.
   */
  void setupDerivativeData();

  /** @brief Initialize all workers and set up derivative metadata. */
  void init();

  /** @brief Compute the LFC for all workers, then update LFC derivatives. */
  void computeLfc();

  /** @brief Compute the SSF for all workers. */
  void computeSsf();

  /**
   * @brief Compute the maximum convergence error across all workers.
   * @return Maximum per-worker RMS difference between successive SSF iterates.
   */
  double computeError() const;

  /** @brief Mix the new and old SSF for all workers. */
  void updateSolution();

  /** @brief Set the initial SSF guess for all workers. */
  void initialGuess();

  /**
   * @brief Return the VS-specific input parameters.
   * @return Constant reference to the @p VSInput.
   */
  virtual const VSInput &inVS() const = 0;

  /**
   * @brief Return the scheme-specific base input parameters.
   * @return Constant reference to the @p Input.
   */
  virtual const Input &inScheme() const = 0;

private:

  /** @brief Compute and store LFC finite-difference derivatives for all grid
   * points. */
  void computeLfcDerivatives();

  /**
   * @brief Compute the finite difference of LFC column @p l at grid index @p i.
   * @param f LFC array.
   * @param l Matsubara frequency column index.
   * @param i Flat grid index.
   * @param t Derivative type (centered / forward / backward).
   * @return Finite-difference value.
   */
  double
  derivative(const Vector2D &f, int l, size_t i, DerivativeData::Type t) const;

  /**
   * @brief Compute a scalar finite difference from three adjacent values.
   * @param f0 "Down" value.
   * @param f1 "Center" value.
   * @param f2 "Up" value.
   * @param t  Derivative type.
   * @return Finite-difference value.
   */
  double
  derivative(double f0, double f1, double f2, DerivativeData::Type t) const;
};

#endif
