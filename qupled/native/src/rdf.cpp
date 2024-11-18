#include "rdf.hpp"

using namespace std;
using ItgParam = Integrator1D::Param;
using ItgType = Integrator1D::Type;

double Rdf::ssf(const double &y) const { return ssfi.eval(y); }

double Rdf::integrand(const double &y) const {
  if (y > cutoff) return 0;
  const double yssf = y * (ssf(y) - 1);
  return (r == 0.0) ? y * yssf : yssf;
}

double Rdf::get() const {
  auto func = [&](const double &y) -> double { return integrand(y); };
  if (r == 0) {
    itg.compute(func, ItgParam(0.0, cutoff));
    return 1 + 1.5 * itg.getSolution();
  } else {
    itgf.compute(func, ItgParam(r));
    return 1 + 1.5 * itgf.getSolution() / r;
  }
}
