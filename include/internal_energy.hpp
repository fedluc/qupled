#ifndef INTERNAL_ENERGY_HPP
#define INTERNAL_ENERGY_HPP

#include "numerics.hpp"
#include <cmath>

// -----------------------------------------------------------------
// Class for internal energy calculation
// -----------------------------------------------------------------

class InternalEnergy {

public:

  // Constructor
  InternalEnergy(const double &rs_,
                 const double &yMin_,
                 const double &yMax_,
                 const Interpolator1D &ssfi_,
                 Integrator1D &itg_)
      : rs(rs_),
        yMin(yMin_),
        yMax(yMax_),
        itg(itg_),
        ssfi(ssfi_) {}

  // Get result of integration
  double get() const;

private:

  // Coupling parameter
  const double rs;

  // Integration limits
  const double yMin;
  const double yMax;

  // Integrator object
  Integrator1D &itg;

  // Static structure factor interpolator
  const Interpolator1D &ssfi;

  // Integrand
  double integrand(const double &y) const;

  // Compute static structure factor
  double ssf(const double &y) const;

  // Constant for unit conversion
  const double lambda = pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
};

#endif
