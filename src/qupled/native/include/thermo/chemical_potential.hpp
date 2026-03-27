#ifndef CHEMICALPOTENTIAL_HPP
#define CHEMICALPOTENTIAL_HPP

#include "schemes/input.hpp"
#include "util/dimensions_util.hpp"
#include <memory>

/**
 * @brief Computes the chemical potential from the density normalization
 * condition.
 *
 * The chemical potential is found by solving the equation
 * @f$\int_0^\infty f(\epsilon; \mu)\, g(\epsilon)\, d\epsilon = n@f$
 * where @f$f@f$ is the Fermi–Dirac distribution and @f$g@f$ is the
 * density of states.  Both 2D and 3D geometries are supported via the
 * @p DimensionsHandler dispatch mechanism.
 */
class ChemicalPotential : public dimensionsUtil::DimensionsHandler {
public:

  /**
   * @brief Construct with the simulation input parameters.
   * @param in_ Shared pointer to the input parameters.
   */
  explicit ChemicalPotential(const std::shared_ptr<const Input> in_)
      : in(in_) {}

  /**
   * @brief Return the computed chemical potential.
   * @return Chemical potential in units of the Fermi energy.
   */
  double get() const { return mu; }

private:

  /** @brief Shared pointer to the input parameters. */
  const std::shared_ptr<const Input> in;

  /** @brief Computed chemical potential (initialized to the sentinel NaN). */
  double mu = DEFAULT_DOUBLE;

  /** @brief Compute the chemical potential for a 2D system. */
  void compute2D() override;

  /** @brief Compute the chemical potential for a 3D system. */
  void compute3D() override;

  /**
   * @brief Evaluate the normalization condition at a trial chemical potential.
   *
   * Used as the root-finding objective function.
   *
   * @param mu Trial chemical potential value.
   * @return Residual: computed density minus target density.
   */
  double normalizationCondition(const double &mu) const;
};

#endif
