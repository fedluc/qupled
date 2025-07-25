#include "rdf.hpp"
#include "gsl/gsl_sf_bessel.h"

using namespace std;
using ItgParam = Integrator1D::Param;
using ItgType = Integrator1D::Type;

double Rdf::ssf(const double &y) const { return ssfi->eval(y); }

double Rdf::integrand(const double &y) const {
  if (y > cutoff) return 0;
  const double yssf = y * (ssf(y) - 1);
  return (r == 0.0) ? y * yssf : yssf;
}

double Rdf::integrand2D(const double &y) const {
  if (y > cutoff) return 0;
  const double yssf = y * (ssf(y) - 1);
  return (r == 0.0) ? yssf : yssf * SpecialFunctions::besselJ0(r * y);
}

double Rdf::get() {
  compute(dim);
  return res;
}

void Rdf::compute2D() {
  auto func = [&](const double &y) -> double { return integrand2D(y); };
  itg->compute(func, ItgParam(0.0, cutoff));
  res = 1 + itg->getSolution();
}

void Rdf::compute3D() {
  auto func = [&](const double &y) -> double { return integrand(y); };
  if (r == 0) {
    itg->compute(func, ItgParam(0.0, cutoff));
    res = 1 + 1.5 * itg->getSolution();
  } else {
    itgf->compute(func, ItgParam(r));
    res = 1 + 1.5 * itgf->getSolution() / r;
  }
}
