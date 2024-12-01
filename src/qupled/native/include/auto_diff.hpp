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
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = dual1.der[i] + dual2.der[i];
  }
  return AutoDiff<T>(dual1.val + dual2.val, der);
}

template <typename T>
AutoDiff<T> operator+(const AutoDiff<T> &dual, const double &scalar) {
  return AutoDiff<T>(dual.val + scalar, dual.der);
}

template <typename T>
AutoDiff<T> operator+(const double &scalar, const AutoDiff<T> &dual) {
  return dual + scalar;
}

// subtraction operators
template <typename T>
AutoDiff<T> operator-(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = dual1.der[i] - dual2.der[i];
  }
  return AutoDiff<T>(dual1.val - dual2.val, der);
}

template <typename T>
AutoDiff<T> operator-(const AutoDiff<T> &dual, double scalar) {
  return AutoDiff<T>(dual.val - scalar, dual.der);
}

template <typename T>
AutoDiff<T> operator-(const double &scalar, const AutoDiff<T> &dual) {
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = -dual.der[i];
  }
  return AutoDiff<T>(scalar - dual.val, der);
}

// multiplication operators
template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = dual1.val * dual2.der[i] + dual1.der[i] * dual2.val;
  }
  return AutoDiff<T>(dual1.val * dual2.val, der);
}

template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual, const double &scalar) {
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = dual.der[i] * scalar;
  }
  return AutoDiff<T>(dual.val * scalar, der);
}

template <typename T>
AutoDiff<T> operator*(const double &scalar, const AutoDiff<T> &dual) {
  return dual * scalar;
}

// division operators
template <typename T>
AutoDiff<T> operator/(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  const double inv_val = 1.0 / dual2.val;
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = (dual1.der[i] - dual1.val * dual2.der[i] * inv_val) * inv_val;
  }
  return AutoDiff<T>(dual1.val * inv_val, der);
}

template <typename T>
AutoDiff<T> operator/(const AutoDiff<T> &dual, const double &scalar) {
  const double inv_scalar = 1.0 / scalar;
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = dual.der[i] * inv_scalar;
  }
  return AutoDiff<T>(dual.val * inv_scalar, der);
}

template <typename T>
AutoDiff<T> operator/(const double &scalar, const AutoDiff<T> &dual) {
  const double inv_val = 1.0 / dual.val;
  const double inv_val2 = inv_val * inv_val;
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = -scalar * dual.der[i] * inv_val2;
  }
  return AutoDiff<T>(scalar * inv_val, der);
}

// exponential function
template <typename T>
AutoDiff<T> exp(const AutoDiff<T> &x) {
  const double exp_val = exp(x.val);
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = exp_val * x.der[i];
  }
  return AutoDiff<T>(exp_val, der);
}

// square root function
template <typename T>
AutoDiff<T> sqrt(const AutoDiff<T> &x) {
  const double sqrt_val = sqrt(x.val);
  const double inv_sqrt = 0.5 / sqrt_val;
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = x.der[i] * inv_sqrt;
  }
  return AutoDiff<T>(sqrt_val, der);
}

// hyperbolic tangent function
template <typename T>
AutoDiff<T> tanh(const AutoDiff<T> &x) {
  const double tanh_val = tanh(x.val);
  const double sech2_val = 1.0 - tanh_val * tanh_val;
  std::vector<T> der(NVAR);
  for (size_t i = 0; i < NVAR; ++i) {
    der[i] = x.der[i] * sech2_val;
  }
  return AutoDiff<T>(tanh_val, der);
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
  std::reference_wrapper<double> dx = der[0];
  std::reference_wrapper<double> dy = der[1];
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
  std::reference_wrapper<double> dx = der[0].val;
  std::reference_wrapper<double> dy = der[1].val;
  std::reference_wrapper<double> dxx = der[0].der[0];
  std::reference_wrapper<double> dxy = der[0].der[1];
  std::reference_wrapper<double> dyy = der[1].der[1];
};

#endif
