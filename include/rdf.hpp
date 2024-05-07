#ifndef RDF_HPP
#define RDF_HPP

#include "numerics.hpp"

// ------------------------------------------------------
// Class for the radial distribution function calculation
// ------------------------------------------------------

class Rdf {

public:

  // Constructor
  Rdf(const double &r_,
      const double &cutoff_,
      const Interpolator1D &ssfi_,
      Integrator1D &itg_,
      Integrator1D &itgf_)
      : r(r_),
        cutoff(cutoff_),
        itgf(itgf_),
        itg(itg_),
        ssfi(ssfi_) {}

  // Get result of integration
  double get() const;

private:

  // Spatial position
  const double r;

  // Cutoff in the wave-vector grid
  const double cutoff;

  // Fourier Integrator object
  Integrator1D &itgf;

  // Integrator object
  Integrator1D &itg;

  // Static structure factor interpolator
  const Interpolator1D &ssfi;

  // Integrand
  double integrand(const double &y) const;

  // Compute static structure factor
  double ssf(const double &y) const;
};

#endif
