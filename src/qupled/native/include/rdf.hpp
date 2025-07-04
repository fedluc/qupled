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
      std::shared_ptr<Interpolator1D> ssfi_,
      std::shared_ptr<Integrator1D> itg_,
      std::shared_ptr<Integrator1D> itgf_)
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
  const std::shared_ptr<Integrator1D> itgf;
  // Integrator object
  const std::shared_ptr<Integrator1D> itg;
  // Static structure factor interpolator
  const std::shared_ptr<Interpolator1D> ssfi;
  // Integrand
  double integrand(const double &y) const;
  // Compute static structure factor
  double ssf(const double &y) const;
};

class Rdf2D {

public:

  // Constructor
  Rdf2D(const double &r_,
      const double &cutoff_,
      std::shared_ptr<Interpolator1D> ssfi_,
      std::shared_ptr<Integrator1D> itg_)
      : r(r_),
        cutoff(cutoff_),
        itg(itg_),
        ssfi(ssfi_) {}

  // Get result of integration
  double get2D() const;

private:

  // Spatial position
  const double r;
  // Cutoff in the wave-vector grid
  const double cutoff;
  // Fourier Integrator object
  const std::shared_ptr<Integrator1D> itgf;
  // Integrator object
  const std::shared_ptr<Integrator1D> itg;
  // Static structure factor interpolator
  const std::shared_ptr<Interpolator1D> ssfi;
  // Integrand
  double integrand2D(const double &y) const;
  // Compute static structure factor
  double ssf2D(const double &y) const;
};

#endif
