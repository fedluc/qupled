#ifndef VECTOR_UTIL_HPP
#define VECTOR_UTIL_HPP

#include <vector>

// ------------------------------------------------------------------
// Utility functions to manipulate vectors from the standard library
// ------------------------------------------------------------------

namespace vecUtil {

  // Element-wise sum between two vectors
  std::vector<double> sum(const std::vector<double> &v1,
                          const std::vector<double> &v2);

  // Element-wise difference between two vectors
  std::vector<double> diff(const std::vector<double> &v1,
                           const std::vector<double> &v2);

  // Element-wise multiplication of two vectors
  std::vector<double> mult(const std::vector<double> &v1,
                           const std::vector<double> &v2);

  // Element-wise division of two vectors
  std::vector<double> div(const std::vector<double> &v1,
                          const std::vector<double> &v2);

  // Element-wise multiplication of a vector and a scalar
  std::vector<double> mult(const std::vector<double> &v, const double a);

  // Linear combination of two vectors
  std::vector<double> linearCombination(const std::vector<double> &v1,
                                        const double a,
                                        const std::vector<double> &v2,
                                        const double b);

  // Root mean square difference between two vectors
  double rms(const std::vector<double> &v1,
             const std::vector<double> &v2,
             const bool normalize);

  // Fill vector with constant values
  void fill(std::vector<double> &v, const double &num);

} // namespace vecUtil

#endif
