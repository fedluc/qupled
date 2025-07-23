#include "internal_energy.hpp"
#include "numerics.hpp"

double InternalEnergy::ssf(const double &y) const { return ssfi->eval(y); }

double InternalEnergy::integrand(const double &y) const { return ssf(y) - 1; }

double InternalEnergy::get() const {
  auto func = [&](const double &y) -> double { return integrand(y); };
  itg->compute(func, Integrator1D::Param(yMin, yMax));
  if (dim == dimensionsUtil::Dimension::D2) {
    return itg->getSolution() * sqrt(2.0) / (2.0 * rs);
  }
  return itg->getSolution() / (M_PI * rs * numUtil::lambda);
}
