#ifndef VS_QVSSTLS_HPP
#define VS_QVSSTLS_HPP

#include "schemes/input.hpp"
#include "schemes/qstls.hpp"
#include "util/numerics.hpp"
#include "vs/vsbase.hpp"
#include "vs/vsmanager.hpp"
#include <memory>

/**
 * @brief qSTLS-based worker for a single grid point in the QVS scheme.
 *
 * Wraps a @p Qstls solver to provide the @p VSWorker interface required by
 * @p VSQstlsManager. Constructed for one (rs, Theta) grid point; all
 * structural property queries delegate to the underlying qSTLS solver.
 */
class VSQstlsWorker : public VSWorker, public Qstls {
public:

  /** @brief Type alias for the input class. */
  using InputType = QVSStlsInput;

  /**
   * @brief Construct a worker for the given input and grid point.
   * @param in Shared pointer to the QVS-STLS input parameters.
   * @param p  Grid point that determines the (rs, Theta) offset for this
   * worker.
   */
  VSQstlsWorker(const std::shared_ptr<const QVSStlsInput> &in, GridPoint p);

  const Vector2D &getLfc() const override { return Qstls::getLfc(); }
  const std::vector<double> &getWvg() const override { return Qstls::getWvg(); }
  const std::vector<double> &getSsf() const override { return Qstls::getSsf(); }
  const Vector2D &getIdr() const override { return Qstls::getIdr(); }
  std::vector<double> getSdr() const override { return Qstls::getSdr(); }
  double getUInt() const override { return Qstls::getUInt(); }

  /**
   * @brief Compute the Q-adder contribution for the quantum VS free-energy
   * expression.
   * @return Q-adder value at this grid point.
   */
  double getQAdder() const override;

  double getCoupling() const override { return inPtr->getCoupling(); }
  double getDegeneracy() const override { return inPtr->getDegeneracy(); }
  void init() override { Qstls::init(); }
  void initialGuess() override { Qstls::initialGuess(); }
  void computeSsf() override { Qstls::computeSsf(); }
  void computeLfc() override { Qstls::computeLfc(); }

  /**
   * @brief Apply a LFC correction by subtracting @p v from the stored LFC.
   * @param v LFC correction to subtract.
   */
  void applyLfcDiff(const Vector2D &v) override { lfc.diff(v); }
  double computeError() const override { return Qstls::computeError(); }
  void updateSolution() override { Qstls::updateSolution(); }

  /**
   * @brief Compute the Q-adder using explicitly provided integrators.
   * @param itg2D   Shared pointer to a 2D numerical integrator.
   * @param itgGrid Grid for 2D integration.
   * @return Q-adder value.
   */
  double computeQAdder(const std::shared_ptr<Integrator2D> &itg2D,
                       const std::vector<double> &itgGrid) const;
};

/**
 * @brief VS manager for qSTLS-based schemes.
 *
 * Combines @p VSManager (grid iteration) with @p Qstls (self-consistency loop)
 * to implement the QVS-STLS variational-swarm scheme. Owns nine
 * @p VSQstlsWorker objects and a shared 2D integrator for the Q-adder
 * computations.
 */
class VSQstlsManager : public VSManager, public Qstls {
public:

  /**
   * @brief Construct the manager and initialize all nine workers.
   * @param in Shared pointer to the QVS-STLS input parameters.
   */
  explicit VSQstlsManager(const std::shared_ptr<const QVSStlsInput> &in);

  void init() override { VSManager::init(); }
  void initialGuess() override { VSManager::initialGuess(); }
  void computeSsf() override { VSManager::computeSsf(); }
  void computeLfc() override { VSManager::computeLfc(); }
  double computeError() const override { return VSManager::computeError(); }
  void updateSolution() override { VSManager::updateSolution(); }

  /**
   * @brief Run the qSTLS self-consistency loop as the manager's compute method.
   * @return 0 on success.
   */
  int compute() override { return Qstls::compute(); }

private:

  /** @brief Shared pointer to the input, used for @p inVS() and @p inScheme().
   */
  std::shared_ptr<const QVSStlsInput> managerInPtr_;
  /** @brief Shared 2D integrator for Q-adder computations. */
  std::shared_ptr<Integrator2D> itg2D;
  /** @brief Grid for 2D integration of the Q-adder. */
  std::vector<double> itgGrid;
  const VSInput &inVS() const override { return *managerInPtr_; }
  const Input &inScheme() const override { return *managerInPtr_; }
};

/**
 * @brief Variational-swarm solver based on qSTLS structural properties.
 *
 * Determines the optimal free parameter @p alpha by minimizing the
 * exchange-correlation free energy over a 3×3 (rs, Theta) grid, using
 * @p VSQstlsManager to run the qSTLS solver at each grid point.
 */
class QVSStls : public VSBase {
public:

  /**
   * @brief Construct the QVS-STLS solver.
   * @param in Shared pointer to the QVS-STLS input parameters.
   */
  explicit QVSStls(const std::shared_ptr<const QVSStlsInput> &in);

  /** @brief Expose @p VSBase::compute() as the public entry point. */
  using VSBase::compute;

private:

  /** @brief Shared pointer to the input parameters. */
  std::shared_ptr<const QVSStlsInput> inPtr;
  /** @brief 3×3 grid manager. */
  VSQstlsManager grid_;

  VSManager &grid() override { return grid_; }
  const VSManager &grid() const override { return grid_; }
  const VSInput &in() const override;
  const Input &inScheme() const override;
};

/**
 * @brief Computes the Q-adder contribution for the quantum VS free-energy
 * expression.
 *
 * The Q-adder is a functional of the static structure factor that enters the
 * self-consistency equation for the free parameter @p alpha in the quantum VS
 * (QVS) scheme. It combines 1D and 2D integrations of the SSF with the
 * Fermi–Dirac kernel.
 */
class QAdder {

public:

  /**
   * @brief Construct the Q-adder calculator.
   * @param Theta_    Degeneracy parameter.
   * @param mu_       Chemical potential.
   * @param limitMin  Lower integration limit.
   * @param limitMax  Upper integration limit.
   * @param itgGrid_  Grid for 2D integration.
   * @param itg1_     Shared pointer to a 1D integrator.
   * @param itg2_     Shared pointer to a 2D integrator.
   * @param interp_   Shared pointer to an SSF interpolator.
   */
  QAdder(const double &Theta_,
         const double &mu_,
         const double &limitMin,
         const double &limitMax,
         const std::vector<double> &itgGrid_,
         std::shared_ptr<Integrator1D> itg1_,
         std::shared_ptr<Integrator2D> itg2_,
         std::shared_ptr<Interpolator1D> interp_)
      : Theta(Theta_),
        mu(mu_),
        limits(limitMin, limitMax),
        itgGrid(itgGrid_),
        itg1(itg1_),
        itg2(itg2_),
        interp(interp_) {}

  /**
   * @brief Compute and return the Q-adder value.
   * @return Q-adder at the current state point.
   */
  double get() const;

private:

  /** @brief Degeneracy parameter. */
  const double Theta;
  /** @brief Chemical potential. */
  const double mu;
  /** @brief Integration limits [min, max]. */
  const std::pair<double, double> limits;
  /** @brief Grid for 2D integration. */
  const std::vector<double> &itgGrid;
  /** @brief 1D numerical integrator. */
  const std::shared_ptr<Integrator1D> itg1;
  /** @brief 2D numerical integrator. */
  const std::shared_ptr<Integrator2D> itg2;
  /** @brief Interpolator for the static structure factor. */
  const std::shared_ptr<Interpolator1D> interp;

  /**
   * @brief Evaluate the interpolated SSF at wave-vector @p y.
   * @param y Wave-vector value.
   */
  double ssf(const double &y) const;

  /**
   * @brief Denominator integrand over @p q.
   * @param q Auxiliary momentum.
   */
  double integrandDenominator(const double q) const;

  /**
   * @brief First numerator integrand over @p q.
   * @param q Auxiliary momentum.
   */
  double integrandNumerator1(const double q) const;

  /**
   * @brief Second numerator integrand over @p w.
   * @param w Auxiliary frequency variable.
   */
  double integrandNumerator2(const double w) const;

  /**
   * @brief Compute the denominator integral.
   * @param res Output variable for the result.
   */
  void getIntDenominator(double &res) const;
};

#endif
