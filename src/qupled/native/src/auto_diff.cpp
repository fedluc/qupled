#include "auto_diff.hpp"
#include <cmath>

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

// addition operators
AutoDiff1 operator+(const AutoDiff1 &dual1, const AutoDiff1 &dual2) {
  return AutoDiff1(
      dual1.val + dual2.val, dual1.dx + dual2.dx, dual1.dy + dual2.dy);
}

AutoDiff1 operator+(const AutoDiff1 &dual, const double &scalar) {
  return AutoDiff1(dual.val + scalar, dual.dx, dual.dy);
}

AutoDiff1 operator+(const double &scalar, const AutoDiff1 &dual) {
  return dual + scalar;
}

// subtraction operators
AutoDiff1 operator-(const AutoDiff1 &dual1, const AutoDiff1 &dual2) {
  return AutoDiff1(
      dual1.val - dual2.val, dual1.dx - dual2.dx, dual1.dy - dual2.dy);
}

AutoDiff1 operator-(const AutoDiff1 &dual, double scalar) {
  return AutoDiff1(dual.val - scalar, dual.dx, dual.dy);
}

AutoDiff1 operator-(const double &scalar, const AutoDiff1 &dual) {
  return AutoDiff1(scalar - dual.val, -dual.dx, -dual.dy);
}

// multiplication operators
AutoDiff1 operator*(const AutoDiff1 &dual1, const AutoDiff1 &dual2) {
  return AutoDiff1(dual1.val * dual2.val,
                   dual1.val * dual2.dx + dual1.dx * dual2.val,
                   dual1.val * dual2.dy + dual1.dy * dual2.val);
}

AutoDiff1 operator*(const AutoDiff1 &dual, const double &scalar) {
  return AutoDiff1(dual.val * scalar, dual.dx * scalar, dual.dy * scalar);
}

AutoDiff1 operator*(const double &scalar, const AutoDiff1 &dual) {
  return dual * scalar;
}

// division operators
AutoDiff1 operator/(const AutoDiff1 &dual1, const AutoDiff1 &dual2) {
  const double inv_val = 1.0 / dual2.val;
  return AutoDiff1(dual1.val * inv_val,
                   (dual1.dx - dual1.val * dual2.dx * inv_val) * inv_val,
                   (dual1.dy - dual1.val * dual2.dy * inv_val) * inv_val);
}

AutoDiff1 operator/(const AutoDiff1 &dual, const double &scalar) {
  const double inv_scalar = 1.0 / scalar;
  return AutoDiff1(
      dual.val * inv_scalar, dual.dx * inv_scalar, dual.dy * inv_scalar);
}

AutoDiff1 operator/(const double &scalar, const AutoDiff1 &dual) {
  const double inv_val = 1.0 / dual.val;
  return AutoDiff1(scalar * inv_val,
                   -scalar * dual.dx * inv_val * inv_val,
                   -scalar * dual.dy * inv_val * inv_val);
}

// exponential function
AutoDiff1 exp(const AutoDiff1 &x) {
  const double exp_val = exp(x.val);
  return AutoDiff1(exp_val, exp_val * x.dx, exp_val * x.dy);
}

// square root function
AutoDiff1 sqrt(const AutoDiff1 &x) {
  const double sqrt_val = sqrt(x.val);
  const double inv_sqrt = 0.5 / sqrt_val;
  return AutoDiff1(sqrt_val, x.dx * inv_sqrt, x.dy * inv_sqrt);
}

// hyperbolic tangent function
AutoDiff1 tanh(const AutoDiff1 &x) {
  const double tanh_val = tanh(x.val);
  const double sech2_val = 1.0 - tanh_val * tanh_val;
  return AutoDiff1(tanh_val, x.dx * sech2_val, x.dy * sech2_val);
}

// addition operators
AutoDiff2 operator+(const AutoDiff2 &dual1, const AutoDiff2 &dual2) {
  return AutoDiff2(
      dual1.val + dual2.val, dual1.dx + dual2.dx, dual1.dy + dual2.dy);
}

AutoDiff2 operator+(const AutoDiff2 &dual, const double &scalar) {
  return AutoDiff2(dual.val + scalar, dual.dx, dual.dy);
}

AutoDiff2 operator+(const double &scalar, const AutoDiff2 &dual) {
  return dual + scalar;
}

// subtraction operators
AutoDiff2 operator-(const AutoDiff2 &dual1, const AutoDiff2 &dual2) {
  return AutoDiff2(
      dual1.val - dual2.val, dual1.dx - dual2.dx, dual1.dy - dual2.dy);
}

AutoDiff2 operator-(const AutoDiff2 &dual, double scalar) {
  return AutoDiff2(dual.val - scalar, dual.dx, dual.dy);
}

AutoDiff2 operator-(const double &scalar, const AutoDiff2 &dual) {
  return AutoDiff2(scalar - dual.val, -1.0 * dual.dx, -1.0 * dual.dy);
}

// multiplication operators
AutoDiff2 operator*(const AutoDiff2 &dual1, const AutoDiff2 &dual2) {
  return AutoDiff2(dual1.val * dual2.val,
                   dual1.val * dual2.dx + dual1.dx * dual2.val,
                   dual1.val * dual2.dy + dual1.dy * dual2.val);
}

AutoDiff2 operator*(const AutoDiff2 &dual, const double &scalar) {
  return AutoDiff2(dual.val * scalar, dual.dx * scalar, dual.dy * scalar);
}

AutoDiff2 operator*(const double &scalar, const AutoDiff2 &dual) {
  return dual * scalar;
}

// division operators
AutoDiff2 operator/(const AutoDiff2 &dual1, const AutoDiff2 &dual2) {
  const auto &inv_val = 1.0 / dual2.val;
  return AutoDiff2(dual1.val * inv_val,
                   (dual1.dx - dual1.val * dual2.dx * inv_val) * inv_val,
                   (dual1.dy - dual1.val * dual2.dy * inv_val) * inv_val);
}

AutoDiff2 operator/(const AutoDiff2 &dual, const double &scalar) {
  const auto &inv_scalar = 1.0 / scalar;
  return AutoDiff2(
      dual.val * inv_scalar, dual.dx * inv_scalar, dual.dy * inv_scalar);
}

AutoDiff2 operator/(const double &scalar, const AutoDiff2 &dual) {
  const auto &inv_val = 1.0 / dual.val;
  return AutoDiff2(scalar * inv_val,
                   -scalar * dual.dx * inv_val * inv_val,
                   -scalar * dual.dy * inv_val * inv_val);
}

// exponential function
AutoDiff2 exp(const AutoDiff2 &x) {
  const auto &exp_val = exp(x.val);
  return AutoDiff2(exp_val, exp_val * x.dx, exp_val * x.dy);
}

// square root function
AutoDiff2 sqrt(const AutoDiff2 &x) {
  const auto &sqrt_val = sqrt(x.val);
  const auto &inv_sqrt = 0.5 / sqrt_val;
  return AutoDiff2(sqrt_val, x.dx * inv_sqrt, x.dy * inv_sqrt);
}

// hyperbolic tangent function
AutoDiff2 tanh(const AutoDiff2 &x) {
  const auto &tanh_val = tanh(x.val);
  const auto &sech2_val = 1.0 - tanh_val * tanh_val;
  return AutoDiff2(tanh_val, x.dx * sech2_val, x.dy * sech2_val);
}
