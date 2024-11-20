#ifndef THERMO_UTIL_HPP
#define THERMO_UTIL_HPP

#include <vector>

// -----------------------------------------------------------------
// Utility functions to compute thermodynamic properties
// -----------------------------------------------------------------

namespace thermoUtil {

  double computeInternalEnergy(const std::vector<double> &wvg,
                               const std::vector<double> &ssf,
                               const double &coupling);

  double computeFreeEnergy(const std::vector<double> &grid,
                           const std::vector<double> &rsu,
                           const double &coupling);

  double computeFreeEnergy(const std::vector<double> &grid,
                           const std::vector<double> &rsu,
                           const double &coupling,
                           const bool normalize);

  std::vector<double> computeRdf(const std::vector<double> &r,
                                 const std::vector<double> &wvg,
                                 const std::vector<double> &ssf);

} // namespace thermoUtil

#endif
