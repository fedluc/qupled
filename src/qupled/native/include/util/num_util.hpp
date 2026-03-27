#ifndef NUM_UTIL_HPP
#define NUM_UTIL_HPP

#include <cmath>
#include <limits>

/**
 * @brief Numeric constants and comparison utilities for double-precision
 * values.
 *
 * Provides sentinel values for uninitialized parameters, a dimensionless
 * constant used throughout the dielectric scheme calculations, and tolerance-
 * based comparison helpers.
 */
namespace numUtil {

  /** @brief Positive infinity sentinel. */
  constexpr double Inf = std::numeric_limits<double>::infinity();

  /** @brief Signaling NaN sentinel for uninitialized double parameters. */
  constexpr double NaN = std::numeric_limits<double>::signaling_NaN();

  /** @brief Signaling NaN cast to int for uninitialized integer parameters. */
  constexpr double iNaN = std::numeric_limits<int>::signaling_NaN();

  /** @brief Tolerance used in floating-point comparisons. */
  constexpr double dtol = 1e-10;

  /** @brief Cube of the Wigner–Seitz constant @f$\lambda^3 = 4/(9\pi)@f$. */
  constexpr double lambda3 = 4.0 / (9.0 * M_PI);

  /** @brief Wigner–Seitz constant @f$\lambda = (4/(9\pi))^{1/3}@f$. */
  const double lambda = pow(lambda3, 1.0 / 3.0);

  /**
   * @brief Test whether @p x is zero within the tolerance @p dtol.
   * @param x Value to test.
   * @return True if @f$|x| \le \text{dtol}@f$.
   */
  bool isZero(const double &x);

  /**
   * @brief Test whether two doubles are equal within the tolerance @p dtol.
   * @param x First value.
   * @param y Second value.
   * @return True if @f$|x - y| \le \text{dtol}@f$.
   */
  bool equalTol(const double &x, const double &y);

  /**
   * @brief Test whether @p x is strictly greater than @p y (within @p dtol).
   * @param x Left-hand operand.
   * @param y Right-hand operand.
   * @return True if @f$x > y + \text{dtol}@f$.
   */
  bool largerThan(const double &x, const double &y);

} // namespace numUtil

#endif
