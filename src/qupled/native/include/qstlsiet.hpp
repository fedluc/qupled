#ifndef QSTLSIET_HPP
#define QSTLSIET_HPP

#include "iet.hpp"
#include "input.hpp"
#include "qstls.hpp"

// -----------------------------------------------------------------
// Solver for the qSTLS-based schemes
// -----------------------------------------------------------------

class QstlsIet : public Qstls {

public:

  // Constructor
  explicit QstlsIet(const std::shared_ptr<const QstlsIetInput> &in_);
  explicit QstlsIet(const QstlsIetInput &in_)
      : QstlsIet(std::make_shared<const QstlsIetInput>(in_)) {};
  // Getters
  const std::vector<double> &getBf() const { return iet.getBf(); }

private:

  // Iet extension
  Iet<QstlsIetInput> iet;
  // Iet contribution to the local field correction
  Vector2D lfcIet;
  // Integrator for 2D integrals
  std::vector<double> itgGrid;
  // Input parameters
  const QstlsIetInput &in() const {
    return *StlsUtil::castInput<Input, QstlsIetInput>(inPtr);
  }
  // Initialize basic properties
  void init() override;
  // Compute auxiliary density response
  void computeLfc() override;
  // Compute auxiliary density response
  void computeAdrFixed();
  // Read initital guess from input
  bool initialGuessFromInput() override;
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
           std::vector<std::shared_ptr<Interpolator1D>> lfci_,
           std::shared_ptr<Interpolator1D> bfi_,
           const std::vector<double> &itgGrid_,
           std::shared_ptr<Integrator2D> itg_)
        : QstlsUtil::AdrBase(Theta_, qMin_, qMax_, x_, ssfi_),
          itg(itg_),
          itgGrid(itgGrid_),
          lfci(lfci_),
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
    const std::vector<std::shared_ptr<Interpolator1D>> lfci;
    // Interpolator for the bridge function contribution
    const std::shared_ptr<Interpolator1D> bfi;
    // Interpolator for the fixed component
    Interpolator2D fixi;
    // Compute dynamic local field correction
    double lfc(const double &y, const int &l) const;
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

} // namespace QstlsIetUtil

#endif
