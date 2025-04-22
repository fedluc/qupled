#ifndef STLSIET_HPP
#define STLSIET_HPP

#include "stls.hpp"

// -----------------------------------------------------------------
// Solver for the STLS scheme
// -----------------------------------------------------------------

class StlsIet : public Stls {

public:

  // Constructors
  explicit StlsIet(const StlsInput &in_)
      : Stls(in_, true) {}
  // Compute scheme
  int compute();
  // Getters
  const std::vector<double> &getBf() const { return bf; }

private:

  // Bridge function (for iet schemes)
  std::vector<double> bf;
  void init();
  // Compute static local field correction
  void computeSlfc();
  // Iterations to solve the stls scheme
  void doIterations();
  // Compute bridge function
  void computeBf();
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

  class BridgeFunction {

  public:

    // Constructor
    BridgeFunction(const std::string &theory_,
                   const std::string &mapping_,
                   const double &rs_,
                   const double &Theta_,
                   const double &x_,
                   std::shared_ptr<Integrator1D> itg_)
        : theory(theory_),
          mapping(mapping_),
          rs(rs_),
          Theta(Theta_),
          x(x_),
          itg(itg_) {}
    // Get result of the integration
    double get() const;

  private:

    // Theory to be solved
    const std::string theory;
    // Iet mapping
    const std::string mapping;
    // Coupling parameter
    const double rs;
    // Degeneracy parameter
    const double Theta;
    // Wave vector
    const double x;
    // Integrator object
    const std::shared_ptr<Integrator1D> itg;
    // Constant for unit conversion
    const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
    // Hypernetted-chain bridge function
    double hnc() const;
    // Ichimaru bridge function
    double ioi() const;
    // Lucco Castello and Tolias bridge function
    double lct() const;
    double lctIntegrand(const double &r, const double &Gamma) const;
    // Coupling parameter to compute the bridge function
    double couplingParameter() const;
  };

} // namespace StlsIetUtil

#endif
