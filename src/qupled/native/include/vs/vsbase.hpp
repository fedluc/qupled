#ifndef VS_VSBASE_HPP
#define VS_VSBASE_HPP

#include "schemes/input.hpp"
#include "util/logger.hpp"
#include "util/numerics.hpp"
#include "vs/grid_point.hpp"
#include <vector>

// Forward declaration
class VSManager;

/**
 * @brief Base class orchestrating the variational-swarm (VS) free-energy
 * optimization.
 *
 * @p VSBase drives the outer loop that determines the free variational
 * parameter
 * @p alpha by minimizing the exchange-correlation free energy. For a given
 * @p alpha it delegates the self-consistency iterations to a @p VSManager,
 * which runs nine (rs, Theta) grid points in parallel. The free-energy
 * integrand is accumulated from the interaction energy at each grid point and
 * integrated over the coupling-parameter grid.
 *
 * Derived classes implement the abstract interface methods (@p in(), @p
 * inScheme(),
 * @p grid()) to supply the scheme-specific input and manager objects.
 */
class VSBase : public Logger {

public:

  /** @brief Default constructor. */
  explicit VSBase() {}

  /** @brief Virtual destructor. */
  virtual ~VSBase() = default;

  /**
   * @brief Run the full VS optimization pipeline.
   * @return 0 on success.
   */
  int compute();

  /**
   * @brief Return the converged value of the free parameter alpha.
   * @return Optimized alpha value.
   */
  double getAlpha() const { return alpha; }

  /**
   * @brief Return the free-energy integrand accumulated during the
   * optimization.
   * @return 2D array indexed by [theta-row][rs-column].
   */
  const std::vector<std::vector<double>> &getFreeEnergyIntegrand() const;

  /**
   * @brief Return the coupling-parameter grid used for the free-energy
   * integration.
   * @return Coupling-parameter grid vector.
   */
  const std::vector<double> &getFreeEnergyGrid() const;

  /**
   * @brief Return the SSF at the central (target) grid point.
   * @return Constant reference to the SSF vector.
   */
  const std::vector<double> &getSsf() const;

  /**
   * @brief Return the ITCF at the central (target) grid point.
   * @return Constant reference to the ITCF array.
   */
  const Vector2D &getItcf() const;

  /**
   * @brief Return the LFC at the central (target) grid point.
   * @return Constant reference to the LFC array.
   */
  const Vector2D &getLfc() const;

  /**
   * @brief Return the wave-vector grid at the central (target) grid point.
   * @return Constant reference to the wave-vector grid.
   */
  const std::vector<double> &getWvg() const;

  /**
   * @brief Return the IDR at the central (target) grid point.
   * @return Constant reference to the IDR array.
   */
  const Vector2D &getIdr() const;

  /**
   * @brief Compute and return the static density response at the central grid
   * point.
   * @return Static density response vector.
   */
  std::vector<double> getSdr() const;

  /**
   * @brief Return the interaction energy at the central grid point.
   * @return Dimensionless interaction energy per particle.
   */
  double getUInt() const;

  /**
   * @brief Return the convergence error from the last iteration.
   * @return RMS difference between successive SSF iterates.
   */
  double getError() const;

protected:

  /** @brief Converged value of the free variational parameter alpha. */
  double alpha;

  /** @brief Coupling-parameter grid for the free-energy integral. */
  std::vector<double> rsGrid;

  /**
   * @brief Free-energy integrand values.
   *
   * Outer index: theta row (DOWN = 0, CENTER = 1, UP = 2).
   * Inner index: rs column (same ordering as @p rsGrid).
   */
  std::vector<std::vector<double>> fxcIntegrand;

  /**
   * @brief Return the VS-specific input parameters.
   * @return Constant reference to @p VSInput.
   */
  virtual const VSInput &in() const = 0;

  /**
   * @brief Return the scheme-specific base input parameters.
   * @return Constant reference to @p Input.
   */
  virtual const Input &inScheme() const = 0;

  /**
   * @brief Return the mutable @p VSManager (non-const).
   * @return Reference to the managed grid.
   */
  virtual VSManager &grid() = 0;

  /**
   * @brief Return the @p VSManager (const).
   * @return Constant reference to the managed grid.
   */
  virtual const VSManager &grid() const = 0;

  /**
   * @brief Run one full self-consistency sweep over all grid points.
   * @return 0 on success.
   */
  int runGrid();

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
   * @brief Return the free-energy integrand value at a given grid point.
   * @param p Grid point identifier.
   * @return Free-energy integrand value.
   */
  double getFxcIntegrandValue(GridPoint p) const;

  /**
   * @brief Return the Q-adder (free-energy derivative) at a given grid point.
   * @param p Grid point identifier.
   * @return Q-adder value.
   */
  double getQAdder(GridPoint p) const;

  /** @brief Build the coupling-parameter grid for the free-energy integration.
   */
  void setRsGrid();

  /** @brief Initialize the free-energy integrand array to the correct shape. */
  void setFxcIntegrand();

  /** @brief Update the free-energy integrand values from the current grid
   * state. */
  void updateFxcIntegrand();

  /**
   * @brief Compute the exchange-correlation free energy at a grid point.
   * @param p         Grid point identifier.
   * @param normalize If true, divide by @f$r_s^2@f$.
   * @return Free-energy value.
   */
  double computeFreeEnergy(GridPoint p, bool normalize) const;

  /**
   * @brief Collect the free-energy integrand column for the target Theta row.
   * @return Vector of integrand values over the rs grid.
   */
  std::vector<double> getFreeEnergyData() const;

  /**
   * @brief Return the grid point corresponding to the target state point.
   * @return @p GridPoints::CENTER.
   */
  GridPoint getOutputGridPoint() const;

private:

  /**
   * @brief Find the optimal alpha by minimizing the free energy.
   * @return Converged alpha value.
   */
  double computeAlpha();

  /**
   * @brief Compute the Q-data vector used in the alpha minimization.
   * @return Q-data values for all rs grid points.
   */
  std::vector<double> computeQData();

  /** @brief Drive the alpha iteration loop to convergence. */
  void doIterations();

  /**
   * @brief Evaluate the residual of the alpha self-consistency equation.
   * @param alphaTmp Trial alpha value.
   * @return Residual value.
   */
  double alphaDifference(const double &alphaTmp);
};

#endif
