#ifndef NUM_UTIL_HPP
#define NUM_UTIL_HPP

#include <cmath>
#include <limits>

// -----------------------------------------------------------------
// Utility functions to handle special cases for double numbers
// -----------------------------------------------------------------

namespace numUtil {

  constexpr double Inf = std::numeric_limits<double>::infinity();
  constexpr double NaN = std::numeric_limits<double>::signaling_NaN();
  constexpr double iNaN = std::numeric_limits<int>::signaling_NaN();
  constexpr double dtol = 1e-10;
  constexpr double lambda3 = 4.0 / (9.0 * M_PI);
  const double lambda = pow(lambda3, 1.0 / 3.0);

  // Check if a double is zero within the tolerance dTol
  bool isZero(const double &x);

  // Compare two doubles within the tolerance in dTol
  bool equalTol(const double &x, const double &y);

  // Check that x > y with a dTol tolerance
  bool largerThan(const double &x, const double &y);

} // namespace numUtil

#endif
