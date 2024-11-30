#ifndef AUTO_DIFF_HPP
#define AUTO_DIFF_HPP

#include <cmath>
#include <vector>

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

template <typename T>
class AutoDiff {};

// First order derivatives
template <>
class AutoDiff<double> {
public:

  double val;
  double dx;
  double dy;

  // Constructors
  AutoDiff(double val_, double dx_ = 0.0, double dy_ = 0.0)
      : val(val_),
        dx(dx_),
        dy(dy_) {}
};
using AutoDiff1 = AutoDiff<double>;

// Second order derivatives
template <>
class AutoDiff<AutoDiff1> {
public:

  AutoDiff1 val;
  AutoDiff1 dx;
  AutoDiff1 dy;

  AutoDiff(AutoDiff1 val_, AutoDiff1 dx_ = 0.0, AutoDiff1 dy_ = 0.0)
      : val(val_),
        dx(dx_),
        dy(dy_) {}

  AutoDiff(const double &val_, const double &dx_ = 0.0, const double &dy_ = 0.0)
      : AutoDiff(AutoDiff1(val_, dx_, dy_), AutoDiff1(dx_), AutoDiff1(dy_)) {}
};
using AutoDiff2 = AutoDiff<AutoDiff1>;

// addition operators
template <typename T>
AutoDiff<T> operator+(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  return AutoDiff<T>(
      dual1.val + dual2.val, dual1.dx + dual2.dx, dual1.dy + dual2.dy);
}

template <typename T>
AutoDiff<T> operator+(const AutoDiff<T> &dual, const double &scalar) {
  return AutoDiff<T>(dual.val + scalar, dual.dx, dual.dy);
}

template <typename T>
AutoDiff<T> operator+(const double &scalar, const AutoDiff<T> &dual) {
  return dual + scalar;
}

// subtraction operators
template <typename T>
AutoDiff<T> operator-(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  return AutoDiff<T>(
      dual1.val - dual2.val, dual1.dx - dual2.dx, dual1.dy - dual2.dy);
}

template <typename T>
AutoDiff<T> operator-(const AutoDiff<T> &dual, double scalar) {
  return AutoDiff<T>(dual.val - scalar, dual.dx, dual.dy);
}

template <typename T>
AutoDiff<T> operator-(const double &scalar, const AutoDiff<T> &dual) {
  return AutoDiff<T>(scalar - dual.val, -dual.dx, -dual.dy);
}

// multiplication operators
template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  return AutoDiff<T>(dual1.val * dual2.val,
                     dual1.val * dual2.dx + dual1.dx * dual2.val,
                     dual1.val * dual2.dy + dual1.dy * dual2.val);
}

template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual, const double &scalar) {
  return AutoDiff<T>(dual.val * scalar, dual.dx * scalar, dual.dy * scalar);
}

template <typename T>
AutoDiff<T> operator*(const double &scalar, const AutoDiff<T> &dual) {
  return dual * scalar;
}

// division operators
template <typename T>
AutoDiff<T> operator/(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  const auto inv_val = 1.0 / dual2.val;
  return AutoDiff<T>(dual1.val * inv_val,
                     (dual1.dx - dual1.val * dual2.dx * inv_val) * inv_val,
                     (dual1.dy - dual1.val * dual2.dy * inv_val) * inv_val);
}

template <typename T>
AutoDiff<T> operator/(const AutoDiff<T> &dual, const double &scalar) {
  const auto inv_scalar = 1.0 / scalar;
  return AutoDiff<T>(
      dual.val * inv_scalar, dual.dx * inv_scalar, dual.dy * inv_scalar);
}

template <typename T>
AutoDiff<T> operator/(const double &scalar, const AutoDiff<T> &dual) {
  const auto inv_val = 1.0 / dual.val;
  return AutoDiff<T>(scalar * inv_val,
                     -scalar * dual.dx * inv_val * inv_val,
                     -scalar * dual.dy * inv_val * inv_val);
}

// exponential function
template <typename T>
AutoDiff<T> exp(const AutoDiff<T> &x) {
  const auto exp_val = exp(x.val);
  return AutoDiff<T>(exp_val, exp_val * x.dx, exp_val * x.dy);
}

// square root function
template <typename T>
AutoDiff<T> sqrt(const AutoDiff<T> &x) {
  const auto sqrt_val = sqrt(x.val);
  const auto inv_sqrt = 0.5 / sqrt_val;
  return AutoDiff<T>(sqrt_val, x.dx * inv_sqrt, x.dy * inv_sqrt);
}

// hyperbolic tangent function
template <typename T>
AutoDiff<T> tanh(const AutoDiff<T> &x) {
  const auto tanh_val = tanh(x.val);
  const auto sech2_val = 1.0 - tanh_val * tanh_val;
  return AutoDiff<T>(tanh_val, x.dx * sech2_val, x.dy * sech2_val);
}

#endif
