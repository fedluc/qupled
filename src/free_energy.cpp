#include "free_energy.hpp"

double FreeEnergy::get() const {
  auto func = [&](const double &y) -> double { return rsui.eval(y); };
  itg.compute(func, Integrator1D::Param(0.0, rs));
  if (normalize) {
    return (rs == 0.0) ? -numUtil::Inf : itg.getSolution() / rs / rs;
  };
  return itg.getSolution();
}
