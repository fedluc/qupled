#ifndef STLSIET_HPP
#define STLSIET_HPP

#include "iet.hpp"
#include "stls.hpp"

// -----------------------------------------------------------------
// Solver for the STLS-IET scheme
// -----------------------------------------------------------------

class StlsIet : public Stls, public Iet {

public:

  // Constructors
  explicit StlsIet(const StlsIetInput &in_)
      : Stls(in_, true),
        Iet(in_, in_, wvg, true),
        in(in_){};

private:

  // Resolve ambiguities
  using Stls::wvg;
  // Input parameters
  StlsIetInput in;
  void init() override;
  // Compute static local field correction
  void computeSlfc() override;
};

namespace StlsIetUtil {

  // -----------------------------------------------------------------
  // Classes for the static local field correction
  // -----------------------------------------------------------------

  class Slfc : public StlsUtil::SlfcBase {

  public:

    // Constructor
    Slfc(const double &x_,
         const double &yMin_,
         const double &yMax_,
         std::shared_ptr<Interpolator1D> ssfi_,
         std::shared_ptr<Interpolator1D> slfci_,
         std::shared_ptr<Interpolator1D> bfi_,
         const std::vector<double> &itgGrid_,
         std::shared_ptr<Integrator2D> itg_)
        : SlfcBase(x_, yMin_, yMax_, ssfi_),
          itg(itg_),
          itgGrid(itgGrid_),
          slfci(slfci_),
          bfi(bfi_) {}
    // Get result of integration
    double get() const;

  private:

    // Integrator object
    const std::shared_ptr<Integrator2D> itg;
    // Grid for 2D integration
    const std::vector<double> itgGrid;
    // Integrands
    double integrand1(const double &y) const;
    double integrand2(const double &w) const;
    // Static local field correction interpolator
    const std::shared_ptr<Interpolator1D> slfci;
    // Bridge function interpolator
    const std::shared_ptr<Interpolator1D> bfi;
    // Compute static local field correction
    double slfc(const double &x) const;
    // Compute bridge function
    double bf(const double &x_) const;
  };

} // namespace StlsIetUtil

#endif
