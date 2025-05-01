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
                 std::shared_ptr<Interpolator1D> ssfi_,
                 std::shared_ptr<Integrator1D> itg_)
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
  const std::shared_ptr<Integrator1D> itg;

  // Static structure factor interpolator
  const std::shared_ptr<Interpolator1D> ssfi;

  // Integrand
  double integrand(const double &y) const;

  // Compute static structure factor
  double ssf(const double &y) const;
};

#endif
