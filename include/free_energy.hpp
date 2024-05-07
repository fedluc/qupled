#ifndef FREE_ENERGY_HPP
#define FREE_ENERGY_HPP

#include "numerics.hpp"

// -----------------------------------------------------------------
// Class for free energy calculation
// -----------------------------------------------------------------

class FreeEnergy {

public:

  // Constructor
  FreeEnergy(const double &rs_,
             const Interpolator1D &rsui_,
             Integrator1D &itg_,
             const bool normalize_)
      : rs(rs_),
        itg(itg_),
        rsui(rsui_),
        normalize(normalize_) {}

  // Get result of integration
  double get() const;

private:

  // Coupling parameter
  const double rs;

  // Integrator object
  Integrator1D &itg;

  // Integrand interpolator (the integrand is given by rs * u)
  const Interpolator1D &rsui;

  // Integrand
  double integrand(const double y) const;

  // Flag marking whether the free energy should be normalized with rs^2
  const bool normalize;
};

#endif
