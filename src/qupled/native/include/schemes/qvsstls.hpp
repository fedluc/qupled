#ifndef VS_QVSSTLS_HPP
#define VS_QVSSTLS_HPP

#include "schemes/input.hpp"
#include "util/numerics.hpp"
#include "schemes/qstls.hpp"
#include "vs/vsbase.hpp"
#include "vs/vsmanager.hpp"
#include <memory>

// -----------------------------------------------------------------
// VSQstlsWorker
// -----------------------------------------------------------------

class VSQstlsWorker : public VSWorker, public Qstls {
public:

  using InputType = QVSStlsInput;

  VSQstlsWorker(const std::shared_ptr<const QVSStlsInput> &in, GridPoint p);

  const Vector2D &getLfc() const override { return Qstls::getLfc(); }
  const std::vector<double> &getWvg() const override { return Qstls::getWvg(); }
  const std::vector<double> &getSsf() const override { return Qstls::getSsf(); }
  const Vector2D &getIdr() const override { return Qstls::getIdr(); }
  std::vector<double> getSdr() const override { return Qstls::getSdr(); }
  double getUInt() const override { return Qstls::getUInt(); }
  double getQAdder() const override;
  double getCoupling() const override { return inPtr->getCoupling(); }
  double getDegeneracy() const override { return inPtr->getDegeneracy(); }
  void init() override { Qstls::init(); }
  void initialGuess() override { Qstls::initialGuess(); }
  void computeSsf() override { Qstls::computeSsf(); }
  void computeLfc() override { Qstls::computeLfc(); }
  void applyLfcDiff(const Vector2D &v) override { lfc.diff(v); }
  double computeError() const override { return Qstls::computeError(); }
  void updateSolution() override { Qstls::updateSolution(); }

  double computeQAdder(const std::shared_ptr<Integrator2D> &itg2D,
                       const std::vector<double> &itgGrid) const;
};

// -----------------------------------------------------------------
// VSQstlsManager
// -----------------------------------------------------------------

class VSQstlsManager : public VSManager, public Qstls {
public:

  explicit VSQstlsManager(const std::shared_ptr<const QVSStlsInput> &in);

  void init() override { VSManager::init(); }
  void initialGuess() override { VSManager::initialGuess(); }
  void computeSsf() override { VSManager::computeSsf(); }
  void computeLfc() override { VSManager::computeLfc(); }
  double computeError() const override { return VSManager::computeError(); }
  void updateSolution() override { VSManager::updateSolution(); }
  int compute() override { return Qstls::compute(); }

private:

  std::shared_ptr<const QVSStlsInput> managerInPtr_;
  std::shared_ptr<Integrator2D> itg2D;
  std::vector<double> itgGrid;
  const VSInput &inVS() const override { return *managerInPtr_; }
  const Input &inScheme() const override { return *managerInPtr_; }
};

// -----------------------------------------------------------------
// QVSStls
// -----------------------------------------------------------------

class QVSStls : public VSBase {
public:

  explicit QVSStls(const std::shared_ptr<const QVSStlsInput> &in);
  using VSBase::compute;

private:

  std::shared_ptr<const QVSStlsInput> inPtr;
  VSQstlsManager grid_;

  VSManager &grid() override { return grid_; }
  const VSManager &grid() const override { return grid_; }
  const VSInput &in() const override;
  const Input &inScheme() const override;
};

// -----------------------------------------------------------------
// Class to handle the Q-adder in the free parameter expression
// -----------------------------------------------------------------

class QAdder {

public:

  // Constructor
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
  // Get Q-adder
  double get() const;

private:

  // Degeneracy parameter
  const double Theta;
  // Chemical potential
  const double mu;
  // Integration limits
  const std::pair<double, double> limits;
  // Grid for 2D integration
  const std::vector<double> &itgGrid;
  // Integrator objects
  const std::shared_ptr<Integrator1D> itg1;
  const std::shared_ptr<Integrator2D> itg2;
  // Interpolator 1D class instance
  const std::shared_ptr<Interpolator1D> interp;

  // SSF interpolation
  double ssf(const double &y) const;
  // Integrands
  double integrandDenominator(const double q) const;
  double integrandNumerator1(const double q) const;
  double integrandNumerator2(const double w) const;
  // Get Integral denominator
  void getIntDenominator(double &res) const;
};
#endif
