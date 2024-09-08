#include "vector_util.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

using namespace std;

namespace vecUtil {

  // Element-wise sum between two vectors
  vector<double> sum(const vector<double> &v1, const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res(v1.size());
    transform(v1.begin(), v1.end(), v2.begin(), res.begin(), plus<double>());
    return res;
  }

  // Element-wise difference between two vectors
  vector<double> diff(const vector<double> &v1, const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res(v1.size());
    transform(v1.begin(), v1.end(), v2.begin(), res.begin(), minus<double>());
    return res;
  }

  // Element-wise multiplication between two vectors
  vector<double> mult(const vector<double> &v1, const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res(v1.size());
    transform(
        v1.begin(), v1.end(), v2.begin(), res.begin(), multiplies<double>());
    return res;
  }

  // Element-wise multiplication between two vectors
  vector<double> div(const vector<double> &v1, const vector<double> &v2) {
    assert(v1.size() == v2.size());
    vector<double> res(v1.size());
    transform(v1.begin(), v1.end(), v2.begin(), res.begin(), divides<double>());
    return res;
  }

  // Element-wise multiplication of a vector and a scalar
  vector<double> mult(const vector<double> &v, const double a) {
    vector<double> res = v;
    transform(
        res.begin(), res.end(), res.begin(), [&a](double c) { return c * a; });
    return res;
  }

  // Element-wise linear combination of two vectors
  std::vector<double> linearCombination(const std::vector<double> &v1,
                                        const double a,
                                        const std::vector<double> &v2,
                                        const double b) {
    const size_t n = v1.size();
    assert(v2.size() == n);
    vector<double> res(n);
    for (size_t i = 0; i < n; ++i) {
      res[i] = a * v1[i] + b * v2[i];
    }
    return res;
  }

  // Root square difference between two vectors
  double rms(const vector<double> &v1,
             const vector<double> &v2,
             const bool normalize) {
    const vector<double> tmp = diff(v1, v2);
    double rms = inner_product(tmp.begin(), tmp.end(), tmp.begin(), 0.0);
    if (normalize) rms /= tmp.size();
    return sqrt(rms);
  }

  // Fill vector with constant value
  void fill(vector<double> &v, const double &num) {
    std::for_each(v.begin(), v.end(), [&](double &vi) { vi = num; });
  }

} // namespace vecUtil
