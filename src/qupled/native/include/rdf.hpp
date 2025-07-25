#ifndef RDF_HPP
#define RDF_HPP

#include "dimensions_util.hpp"
#include "input.hpp"
#include "numerics.hpp"

// ------------------------------------------------------
// Class for the radial distribution function calculation
// ------------------------------------------------------

class Rdf : public dimensionsUtil::DimensionsHandler {

public:

  // Constructor
  Rdf(const double &r_,
      const double &cutoff_,
      std::shared_ptr<Interpolator1D> ssfi_,
      std::shared_ptr<Integrator1D> itg_,
      std::shared_ptr<Integrator1D> itgf_,
      const dimensionsUtil::Dimension &dim_)
      : r(r_),
        cutoff(cutoff_),
        itgf(itgf_),
        itg(itg_),
        ssfi(ssfi_),
        dim(dim_),
        res(numUtil::NaN) {}

  // Get result of integration
  double get();

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
  // Dimension of the system
  const dimensionsUtil::Dimension dim;
  // Result of integration
  double res;
  // Integrand
  double integrand(const double &y) const;
  double integrand2D(const double &y) const;
  // Compute static structure factor
  double ssf(const double &y) const;
  // Compute methods
  void compute2D() override;
  void compute3D() override;
};

#endif
