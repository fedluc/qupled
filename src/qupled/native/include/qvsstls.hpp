#ifndef QVSSTLS_HPP
#define QVSSTLS_HPP

#include "input.hpp"
#include "numerics.hpp"
#include "qstls.hpp"
#include "state_point_grid.hpp"
#include "vsbase.hpp"
#include <memory>

// -----------------------------------------------------------------
// VSQstlsWorker: Qstls extended with VS-specific methods for StatePointGrid
// -----------------------------------------------------------------

class VSQstlsWorker : public Qstls {

public:

  using InputType = QVSStlsInput;

  VSQstlsWorker(const std::shared_ptr<const QVSStlsInput> &in,
                bool verbose,
                GridPoint p);

  // VS-specific: apply alpha-correction to lfc (called by StatePointGrid step 3)
  void applyLfcDiff(const Vector2D &v) { lfc.diff(v); }

  // Expose protected Qstls/Stls iteration steps for StatePointGrid orchestration
  void   doInit()           { Qstls::init(); }
  void   doInitialGuess()   { Stls::initialGuess(); }
  void   doComputeLfc()     { Qstls::computeLfc(); }
  void   doComputeSsf()     { Stls::computeSsf(); }
  void   doUpdateSolution() { Stls::updateSolution(); }
  double doComputeError() const { return Stls::computeError(); }

  // Compute quantum Q-adder value using internal mu, wvg, ssf
  double computeQAdder(const std::shared_ptr<Integrator2D> &itg2D,
                       const std::vector<double> &itgGrid) const;
};

// -----------------------------------------------------------------
// QAdder: computes the quantum Q-term in the VS free parameter expression
// -----------------------------------------------------------------

class QAdder {

public:

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

  double get() const;

private:

  const double                      Theta;
  const double                      mu;
  const std::pair<double, double>   limits;
  const std::vector<double>        &itgGrid;
  const std::shared_ptr<Integrator1D>   itg1;
  const std::shared_ptr<Integrator2D>   itg2;
  const std::shared_ptr<Interpolator1D> interp;

  double ssf(const double &y) const;
  double integrandDenominator(const double q) const;
  double integrandNumerator1(const double q) const;
  double integrandNumerator2(const double w) const;
  void   getIntDenominator(double &res) const;
};

// -----------------------------------------------------------------
// QVSStls: VS scheme based on Qstls structural properties
// -----------------------------------------------------------------

class QVSStls : public VSBase {

public:

  explicit QVSStls(const std::shared_ptr<const QVSStlsInput> &in);
  using VSBase::compute;

  // Public interface required by Python bindings
  const std::vector<double> &getSsf()  const { return ssf; }
  const Vector2D &            getLfc()  const { return lfc; }
  const std::vector<double> &getWvg()  const;
  const Vector2D &            getIdr()  const;
  std::vector<double>         getSdr()  const;
  double                      getUInt() const;
  double                      getError() const;

private:

  std::shared_ptr<const QVSStlsInput> inPtr;
  StatePointGrid<VSQstlsWorker>       grid;
  std::shared_ptr<Integrator2D>       itg2D;
  std::vector<double>                 itgGrid;
  std::vector<double>                 ssf;
  Vector2D                            lfc;

  const VSInput &in()       const override;
  const Input   &inScheme() const override;
  void init()           override;
  void updateSolution() override;

  int    runGrid()                               override;
  double getCoupling(GridPoint p)    const override;
  double getDegeneracy(GridPoint p)  const override;
  double getFxcIntegrandValue(GridPoint p) const override;
  // Returns {Q, Qr, Qt} — quantum Q-term and derivatives
  std::vector<double> computeQData() override;
};

#endif
