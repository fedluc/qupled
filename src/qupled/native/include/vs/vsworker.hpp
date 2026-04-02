#ifndef VS_VSWORKER_HPP
#define VS_VSWORKER_HPP

#include "util/vector2D.hpp"
#include <vector>

/**
 * @brief Abstract interface for a single-point variational-swarm worker.
 *
 * Each @p VSWorker computes structural and thermodynamic properties for one
 * (rs, Theta) point in the variational-swarm 3×3 grid. @p VSManager holds
 * an array of nine workers and coordinates the synchronized local-field-
 * correction update that mixes derivatives across grid points.
 *
 * Concrete workers are created by the scheme-specific @p VSManager subclass
 * (e.g., one wrapping a @p Stls solver, another wrapping a @p Qstls solver).
 */
class VSWorker {
public:

  /** @brief Virtual destructor. */
  virtual ~VSWorker() = default;

  /** @brief Return the local field correction (wave-vectors × Matsubara
   * frequencies). */
  virtual const Vector2D &getLfc() const = 0;

  /** @brief Return the wave-vector grid. */
  virtual const std::vector<double> &getWvg() const = 0;

  /** @brief Return the static structure factor over the wave-vector grid. */
  virtual const std::vector<double> &getSsf() const = 0;

  /** @brief Return the ideal density response. */
  virtual const Vector2D &getIdr() const = 0;

  /**
   * @brief Compute and return the static density response.
   * @return Vector of static density response values.
   */
  virtual std::vector<double> getSdr() const = 0;

  /**
   * @brief Return the interaction energy per particle.
   * @return Dimensionless interaction energy.
   */
  virtual double getUInt() const = 0;

  /**
   * @brief Return the chemical potential.
   * @return Chemical potential (in units of the thermal energy).
   */
  virtual double getChemicalPotential() const = 0;

  /**
   * @brief Return the free-energy integrand contribution (Q-adder).
   * @return Q-adder value at this grid point.
   */
  virtual double getQAdder() const = 0;

  /**
   * @brief Return the coupling parameter for this grid point.
   * @return rs value.
   */
  virtual double getCoupling() const = 0;

  /**
   * @brief Return the degeneracy parameter for this grid point.
   * @return Theta value.
   */
  virtual double getDegeneracy() const = 0;

  /** @brief Initialize the worker (wave-vector grid, chemical potential, etc.).
   */
  virtual void init() = 0;

  /** @brief Set the initial guess for the iterative solver. */
  virtual void initialGuess() = 0;

  /** @brief Compute the static structure factor using the current LFC. */
  virtual void computeSsf() = 0;

  /** @brief Compute the local field correction using the current SSF. */
  virtual void computeLfc() = 0;

  /**
   * @brief Apply a corrected LFC difference to update the worker's LFC.
   *
   * Used by @p VSManager to propagate the finite-difference-corrected LFC
   * back to each worker after the synchronized derivative computation.
   *
   * @param v LFC correction array to add.
   */
  virtual void applyLfcDiff(const Vector2D &v) = 0;

  /**
   * @brief Compute the convergence error between successive SSF iterates.
   * @return RMS difference between the current and previous SSF.
   */
  virtual double computeError() const = 0;

  /** @brief Mix the new and old SSF to advance the iteration. */
  virtual void updateSolution() = 0;
};

#endif
