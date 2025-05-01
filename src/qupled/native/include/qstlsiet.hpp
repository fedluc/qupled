#ifndef QSTLSIET_HPP
#define QSTLSIET_HPP

#include "iet.hpp"
#include "input.hpp"
#include "qstls.hpp"

// -----------------------------------------------------------------
// Solver for the qSTLS-based schemes
// -----------------------------------------------------------------

class QstlsIet : public Qstls, public Iet {

public:

  // Constructor
  explicit QstlsIet(const QstlsIetInput &in_);
  // Compute qstls scheme
  int compute();

private:

  // Resolve ambiguities
  using Qstls::wvg;
  // Input data
  const QstlsIetInput in;
  // Auxiliary density response
  Vector2D adrOld;
  // Initialize basic properties
  void init();
  // Compute auxiliary density response
  void computeAdr();
  // Compute static structure factor at finite temperature
  void computeSsf();
  // Iterations to solve the stls scheme
  void doIterations();
  void initialGuess();
  void updateSolution();
  // Compute auxiliary density response
  void computeAdrFixed();
  // Iterations to solve the stls scheme
  bool initialGuessFromInput();
};

namespace QstlsIetUtil {

  class AdrIet : public QstlsUtil::AdrBase {

  public:

    // Constructor for finite temperature calculations
    AdrIet(const double &Theta_,
           const double &qMin_,
           const double &qMax_,
           const double &x_,
           std::shared_ptr<Interpolator1D> ssfi_,
           std::vector<std::shared_ptr<Interpolator1D>> dlfci_,
           std::shared_ptr<Interpolator1D> bfi_,
           const std::vector<double> &itgGrid_,
           std::shared_ptr<Integrator2D> itg_)
        : QstlsUtil::AdrBase(Theta_, qMin_, qMax_, x_, ssfi_),
          itg(itg_),
          itgGrid(itgGrid_),
          dlfci(dlfci_),
          bfi(bfi_) {}

    // Get integration result
    void
    get(const std::vector<double> &wvg, const Vector3D &fixed, Vector2D &res);

  private:

    // Integration limits
    const double &qMin = yMin;
    const double &qMax = yMax;
    // Integrands
    double integrand1(const double &q, const int &l) const;
    double integrand2(const double &y) const;
    // Integrator object
    const std::shared_ptr<Integrator2D> itg;
    // Grid for 2D integration
    const std::vector<double> &itgGrid;
    // Interpolator for the dynamic local field correction
    const std::vector<std::shared_ptr<Interpolator1D>> dlfci;
    // Interpolator for the bridge function contribution
    const std::shared_ptr<Interpolator1D> bfi;
    // Interpolator for the fixed component
    Interpolator2D fixi;
    // Compute dynamic local field correction
    double dlfc(const double &y, const int &l) const;
    // Compute bridge function contribution
    double bf(const double &y) const;
    // Compute fixed component
    double fix(const double &x, const double &y) const;
  };

  class AdrFixedIet : public QstlsUtil::AdrFixedBase {

  public:

    // Constructor for finite temperature calculations
    AdrFixedIet(const double &Theta_,
                const double &qMin_,
                const double &qMax_,
                const double &x_,
                const double &mu_,
                std::shared_ptr<Integrator1D> itg_)
        : QstlsUtil::AdrFixedBase(Theta_, qMin_, qMax_, x_, mu_),
          itg(itg_) {}

    // Get integration result
    void get(int l, const std::vector<double> &wvg, Vector3D &res) const;

  private:

    // Integration limits
    const double &tMin = qMin;
    const double &tMax = qMax;
    // Integrands
    double integrand(const double &t,
                     const double &y,
                     const double &q,
                     const double &l) const;
    // Integrator object
    const std::shared_ptr<Integrator1D> itg;
  };

  // -----------------------------------------------------------------
  // Class for the static structure factor
  // -----------------------------------------------------------------

  class Ssf : public QstlsUtil::Ssf {

  public:

    // Constructor for quantum schemes
    Ssf(const double &x_,
        const double &Theta_,
        const double &rs_,
        const double &ssfHF_,
        std::span<const double> idr_,
        std::span<const double> adr_,
        const double &bf_)
        : QstlsUtil::Ssf(x_, Theta_, rs_, ssfHF_, idr_, adr_),
          bf(bf_) {}
    // Get static structure factor
    double get() const;

  private:

    // Bridge function
    const double bf;
  };

} // namespace QstlsIetUtil

#endif
