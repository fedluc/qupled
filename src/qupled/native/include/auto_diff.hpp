#ifndef AUTO_DIFF_HPP
#define AUTO_DIFF_HPP

#include <cmath>
#include <vector>

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

constexpr size_t NVAR = 2;

template <typename T>
class AutoDiff {
public:

  double val;
  std::vector<T> der;
  // Constructors
  AutoDiff(const double &val_, const std::vector<T> &der_)
      : val(val_),
        der(der_) {}
  explicit AutoDiff(const double &val_)
      : val(val_),
        der(NVAR) {}
};

// addition operators
template <typename T>
AutoDiff<T> operator+(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  AutoDiff<T> result = AutoDiff<T>(dual1.val + dual2.val);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = dual1.der[i] + dual2.der[i];
  }
  return result;
}

template <typename T>
AutoDiff<T> operator+(const AutoDiff<T> &dual, const double &scalar) {
  AutoDiff<T> result = AutoDiff<T>(dual.val + scalar);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = dual.der[i];
  }
  return result;
}

template <typename T>
AutoDiff<T> operator+(const double &scalar, const AutoDiff<T> &dual) {
  return dual + scalar;
}

// subtraction operators
template <typename T>
AutoDiff<T> operator-(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  AutoDiff<T> result = AutoDiff<T>(dual1.val - dual2.val);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = dual1.der[i] - dual2.der[i];
  }
  return result;
}

template <typename T>
AutoDiff<T> operator-(const AutoDiff<T> &dual, double scalar) {
  AutoDiff<T> result = AutoDiff<T>(dual.val - scalar);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = dual.der[i];
  }
  return result;
}

template <typename T>
AutoDiff<T> operator-(const double &scalar, const AutoDiff<T> &dual) {
  AutoDiff<T> result = AutoDiff<T>(scalar - dual.val);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = -dual.der[i];
  }
  return result;
}

// multiplication operators
template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  AutoDiff<T> result = AutoDiff<T>(dual1.val * dual2.val);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = dual1.val * dual2.der[i] + dual1.der[i] * dual2.val;
  }
  return result;
}

template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual, const double &scalar) {
  AutoDiff<T> result = AutoDiff<T>(dual.val * scalar);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = dual.der[i] * scalar;
  }
  return result;
}

template <typename T>
AutoDiff<T> operator*(const double &scalar, const AutoDiff<T> &dual) {
  return dual * scalar;
}

// division operators
template <typename T>
AutoDiff<T> operator/(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  const double inv_val = 1.0 / dual2.val;
  AutoDiff<T> result = AutoDiff<T>(dual1.val * inv_val);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] =
        (dual1.der[i] - dual1.val * dual2.der[i] * inv_val) * inv_val;
  }
  return result;
}

template <typename T>
AutoDiff<T> operator/(const AutoDiff<T> &dual, const double &scalar) {
  const double inv_scalar = 1.0 / scalar;
  AutoDiff<T> result = AutoDiff<T>(dual.val * inv_scalar);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = dual.der[i] * inv_scalar;
  }
  return result;
}

template <typename T>
AutoDiff<T> operator/(const double &scalar, const AutoDiff<T> &dual) {
  const double inv_val = 1.0 / dual.val;
  const double inv_val2 = inv_val * inv_val;
  AutoDiff<T> result = AutoDiff<T>(scalar * inv_val);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = -scalar * dual.der[i] * inv_val2;
  }
  return result;
}

// exponential function
template <typename T>
AutoDiff<T> exp(const AutoDiff<T> &x) {
  const double exp_val = exp(x.val);
  AutoDiff<T> result = AutoDiff<T>(exp_val);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = exp_val * x.der[i];
  }
  return result;
}

// square root function
template <typename T>
AutoDiff<T> sqrt(const AutoDiff<T> &x) {
  const double sqrt_val = sqrt(x.val);
  const double inv_sqrt = 0.5 / sqrt_val;
  AutoDiff<T> result = AutoDiff<T>(sqrt_val);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = x.der[i] * inv_sqrt;
  }
  return result;
}

// hyperbolic tangent function
template <typename T>
AutoDiff<T> tanh(const AutoDiff<T> &x) {
  const double tanh_val = tanh(x.val);
  const double sech2_val = 1.0 - tanh_val * tanh_val;
  AutoDiff<T> result = AutoDiff<T>(tanh_val);
  for (size_t i = 0; i < NVAR; ++i) {
    result.der[i] = x.der[i] * sech2_val;
  }
  return result;
}

// First order derivatives
class AutoDiff1 : public AutoDiff<double> {
public:

  // Constructors
  AutoDiff1(double val_, std::vector<double> der_)
      : AutoDiff<double>(val_, der_) {}
  explicit AutoDiff1(double val_ = 0.0)
      : AutoDiff<double>(val_, std::vector<double>(2, 0)) {}
  AutoDiff1(const AutoDiff<double> &other)
      : AutoDiff<double>(other) {}
  // Aliases for convenient access of the results
  double dx = der[0];
  double dy = der[1];
};

// Second order derivatives
class AutoDiff2 : public AutoDiff<AutoDiff1> {
public:

  // Constructors
  AutoDiff2(double val_, std::vector<AutoDiff1> der_)
      : AutoDiff<AutoDiff1>(val_, der_) {}
  explicit AutoDiff2(double val_)
      : AutoDiff<AutoDiff1>(val_, std::vector<AutoDiff1>(2, AutoDiff1(0.0))) {}
  AutoDiff2(const AutoDiff<AutoDiff1> &other)
      : AutoDiff<AutoDiff1>(other) {}
  // Aliases for convenient access of the results
  double dx = der[0].val;
  double dy = der[1].val;
  double dxx = der[0].der[0];
  double dxy = der[0].der[1];
  double dyy = der[1].der[1];
};

#endif
