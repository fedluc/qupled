#ifndef RDF_HPP
#define RDF_HPP

#include "schemes/input.hpp"
#include "util/dimensions_util.hpp"
#include "util/numerics.hpp"

/**
 * @brief Computes the radial distribution function (RDF) at a single real-space
 * distance.
 *
 * The RDF is obtained by a Fourier transform of the static structure factor:
 * @f$g(r) = 1 + \frac{1}{n} \int \frac{d^d k}{(2\pi)^d}\,
 * e^{i\mathbf{k}\cdot\mathbf{r}} [S(k) - 1]@f$. Both 2D and 3D geometries are
 * supported via the @p DimensionsHandler dispatch.
 */
class Rdf : public dimensionsUtil::DimensionsHandler {

public:

  /**
   * @brief Construct for evaluation at a single real-space point.
   * @param r_      Real-space distance at which to evaluate g(r).
   * @param cutoff_ Wave-vector cutoff (upper limit of the SSF grid).
   * @param ssfi_   Shared pointer to an SSF interpolator.
   * @param itg_    Shared pointer to a 1D numerical integrator (standard
   * integrals).
   * @param itgf_   Shared pointer to a 1D Fourier integrator.
   * @param dim_    Spatial dimension of the system.
   */
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

  /**
   * @brief Compute and return g(r) at the configured distance.
   * @return Radial distribution function value at @p r.
   */
  double get();

private:

  /** @brief Real-space distance. */
  const double r;
  /** @brief Wave-vector cutoff (upper integration limit). */
  const double cutoff;
  /** @brief Fourier-type 1D integrator. */
  const std::shared_ptr<Integrator1D> itgf;
  /** @brief Standard 1D integrator. */
  const std::shared_ptr<Integrator1D> itg;
  /** @brief Interpolator for the static structure factor. */
  const std::shared_ptr<Interpolator1D> ssfi;
  /** @brief Spatial dimension of the system. */
  const dimensionsUtil::Dimension dim;
  /** @brief Result of the RDF computation. */
  double res;

  /**
   * @brief 3D integrand @f$k^2 [S(k) - 1] \sin(kr) / (kr)@f$.
   * @param y Wave-vector value.
   */
  double integrand(const double &y) const;

  /**
   * @brief 2D integrand @f$k [S(k) - 1] J_0(kr)@f$.
   * @param y Wave-vector value.
   */
  double integrand2D(const double &y) const;

  /**
   * @brief Evaluate the interpolated SSF at wave-vector @p y.
   * @param y Wave-vector value.
   * @return Interpolated SSF value.
   */
  double ssf(const double &y) const;

  /** @brief Compute the RDF for a 2D system. */
  void compute2D() override;

  /** @brief Compute the RDF for a 3D system. */
  void compute3D() override;
};

#endif
