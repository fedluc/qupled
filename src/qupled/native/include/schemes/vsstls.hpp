#ifndef VS_VSSTLS_HPP
#define VS_VSSTLS_HPP

#include "schemes/input.hpp"
#include "schemes/stls.hpp"
#include "thermo/thermo_util.hpp"
#include "vs/vsbase.hpp"
#include "vs/vsmanager.hpp"

/**
 * @brief STLS-based worker for a single grid point in the VS scheme.
 *
 * Wraps a @p Stls solver to provide the @p VSWorker interface required by
 * @p VSManager. Each instance is constructed for one (rs, Theta) point of
 * the 3×3 variational-swarm grid and delegates all structural property
 * computations to the underlying STLS solver.
 */
class VSStlsWorker : public VSWorker, public Stls {
public:

  /**
   * @brief Construct a worker for the given input parameters.
   * @param in Shared pointer to the VS-STLS input (verbosity disabled).
   */
  VSStlsWorker(const std::shared_ptr<const VSStlsInput> &in)
      : Stls(in, false) {}

  const Vector2D &getLfc() const override { return Stls::getLfc(); }
  const std::vector<double> &getWvg() const override { return Stls::getWvg(); }
  const std::vector<double> &getSsf() const override { return Stls::getSsf(); }
  const Vector2D &getItcf() const override { return Stls::getItcf(); }
  const Vector2D &getIdr() const override { return Stls::getIdr(); }
  std::vector<double> getSdr() const override { return Stls::getSdr(); }
  double getUInt() const override { return Stls::getUInt(); }

  /**
   * @brief Compute the Q-adder contribution (interaction energy at rs = 1).
   * @return Internal energy per particle evaluated at unit coupling strength.
   */
  double getQAdder() const override {
    return thermoUtil::computeInternalEnergy(
        wvg, ssf, 1.0, inPtr->getDimension());
  }

  double getCoupling() const override { return inPtr->getCoupling(); }
  double getDegeneracy() const override { return inPtr->getDegeneracy(); }
  void init() override { Stls::init(); }
  void initialGuess() override { Stls::initialGuess(); }
  void computeSsf() override { Stls::computeSsf(); }
  void computeLfc() override { Stls::computeLfc(); }

  /**
   * @brief Apply a LFC correction by subtracting @p v from the stored LFC.
   * @param v LFC correction to subtract.
   */
  void applyLfcDiff(const Vector2D &v) override { lfc.diff(v); }
  double computeError() const override { return Stls::computeError(); }
  void updateSolution() override { Stls::updateSolution(); }
};

/**
 * @brief VS manager for STLS-based schemes.
 *
 * Combines @p VSManager (which drives the 3×3 grid iteration) with @p Stls
 * (whose @p compute() runs the self-consistency loop for the central grid
 * point). The nine @p VSStlsWorker objects are initialized in the constructor
 * and coordinated through the @p VSManager interface.
 */
class VSStlsManager : public VSManager, public Stls {
public:

  /**
   * @brief Construct the manager and initialize all nine workers.
   * @param in Shared pointer to the VS-STLS input parameters.
   */
  explicit VSStlsManager(const std::shared_ptr<const VSStlsInput> &in);

  void init() override { VSManager::init(); }
  void initialGuess() override { VSManager::initialGuess(); }
  void computeSsf() override { VSManager::computeSsf(); }
  void computeLfc() override { VSManager::computeLfc(); }
  double computeError() const override { return VSManager::computeError(); }
  void updateSolution() override { VSManager::updateSolution(); }

  /**
   * @brief Run the STLS self-consistency loop as the manager's compute method.
   * @return 0 on success.
   */
  int compute() override { return Stls::compute(); }

private:

  /** @brief Shared pointer to the input, used for @p inVS() and @p inScheme().
   */
  std::shared_ptr<const VSStlsInput> managerInPtr;
  const VSInput &inVS() const override { return *managerInPtr; }
  const Input &inScheme() const override { return *managerInPtr; }
};

/**
 * @brief Variational-swarm solver based on STLS structural properties.
 *
 * Determines the optimal free parameter @p alpha by minimizing the
 * exchange-correlation free energy over a 3×3 (rs, Theta) grid, using
 * @p VSStlsManager to run the STLS solver at each grid point.
 */
class VSStls : public VSBase {
public:

  /**
   * @brief Construct the VS-STLS solver.
   * @param in Shared pointer to the VS-STLS input parameters.
   */
  explicit VSStls(const std::shared_ptr<const VSStlsInput> &in);

  /** @brief Expose @p VSBase::compute() as the public entry point. */
  using VSBase::compute;

private:

  /** @brief Shared pointer to the input parameters. */
  std::shared_ptr<const VSStlsInput> inPtr;
  /** @brief 3×3 grid manager. */
  VSStlsManager grid_;

  VSManager &grid() override { return grid_; }
  const VSManager &grid() const override { return grid_; }
  const VSInput &in() const override;
  const Input &inScheme() const override;
};

#endif
