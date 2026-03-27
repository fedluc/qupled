#ifndef VECTOR_UTIL_HPP
#define VECTOR_UTIL_HPP

#include <vector>

/**
 * @brief Free-function utilities for element-wise operations on @c
 * std::vector<double>.
 *
 * All binary functions require the two input vectors to have the same length
 * and return a new vector of that length.
 */
namespace vecUtil {

  /**
   * @brief Compute the element-wise sum of two vectors.
   * @param v1 First operand.
   * @param v2 Second operand (same length as @p v1).
   * @return Vector @f$v_1 + v_2@f$.
   */
  std::vector<double> sum(const std::vector<double> &v1,
                          const std::vector<double> &v2);

  /**
   * @brief Compute the element-wise difference @f$v_1 - v_2@f$.
   * @param v1 Minuend.
   * @param v2 Subtrahend (same length as @p v1).
   * @return Vector @f$v_1 - v_2@f$.
   */
  std::vector<double> diff(const std::vector<double> &v1,
                           const std::vector<double> &v2);

  /**
   * @brief Compute the element-wise product @f$v_1 \cdot v_2@f$.
   * @param v1 First operand.
   * @param v2 Second operand (same length as @p v1).
   * @return Vector @f$v_1 \odot v_2@f$.
   */
  std::vector<double> mult(const std::vector<double> &v1,
                           const std::vector<double> &v2);

  /**
   * @brief Compute the element-wise quotient @f$v_1 / v_2@f$.
   * @param v1 Dividend.
   * @param v2 Divisor (same length as @p v1; no zero-division check).
   * @return Vector @f$v_1 \oslash v_2@f$.
   */
  std::vector<double> div(const std::vector<double> &v1,
                          const std::vector<double> &v2);

  /**
   * @brief Scale a vector by a scalar.
   * @param v Vector to scale.
   * @param a Scalar factor.
   * @return Vector @f$a \cdot v@f$.
   */
  std::vector<double> mult(const std::vector<double> &v, const double a);

  /**
   * @brief Compute the linear combination @f$a \cdot v_1 + b \cdot v_2@f$.
   * @param v1 First vector.
   * @param a  Coefficient for @p v1.
   * @param v2 Second vector (same length as @p v1).
   * @param b  Coefficient for @p v2.
   * @return Linear combination vector.
   */
  std::vector<double> linearCombination(const std::vector<double> &v1,
                                        const double a,
                                        const std::vector<double> &v2,
                                        const double b);

  /**
   * @brief Compute the root-mean-square difference between two vectors.
   * @param v1        First vector.
   * @param v2        Second vector (same length as @p v1).
   * @param normalize If true, divide the RMS by the RMS of @p v1 (relative
   * error).
   * @return RMS (or relative RMS) difference.
   */
  double rms(const std::vector<double> &v1,
             const std::vector<double> &v2,
             const bool normalize);

  /**
   * @brief Fill every element of @p v with the constant @p num.
   * @param v   Vector to fill (modified in place).
   * @param num Fill value.
   */
  void fill(std::vector<double> &v, const double &num);

} // namespace vecUtil

#endif
