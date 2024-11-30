#ifndef AUTO_DIFF_HPP
#define AUTO_DIFF_HPP

#include <cmath>
#include <vector>

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

template <typename T>
class AutoDiff {
public:

  T val_;
  T dx_;
  T dy_;

  // Constructors
  AutoDiff(const T &val__, const T &dx__, const T &dy__)
      : val_(val__),
        dx_(dx__),
        dy_(dy__) {}
};

// addition operators
template <typename T>
AutoDiff<T> operator+(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  return AutoDiff<T>(
      dual1.val_ + dual2.val_, dual1.dx_ + dual2.dx_, dual1.dy_ + dual2.dy_);
}

template <typename T>
AutoDiff<T> operator+(const AutoDiff<T> &dual, const double &scalar) {
  return AutoDiff<T>(dual.val_ + scalar, dual.dx_, dual.dy_);
}

template <typename T>
AutoDiff<T> operator+(const double &scalar, const AutoDiff<T> &dual) {
  return dual + scalar;
}

// subtraction operators
template <typename T>
AutoDiff<T> operator-(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  return AutoDiff<T>(
      dual1.val_ - dual2.val_, dual1.dx_ - dual2.dx_, dual1.dy_ - dual2.dy_);
}

template <typename T>
AutoDiff<T> operator-(const AutoDiff<T> &dual, double scalar) {
  return AutoDiff<T>(dual.val_ - scalar, dual.dx_, dual.dy_);
}

template <typename T>
AutoDiff<T> operator-(const double &scalar, const AutoDiff<T> &dual) {
  return AutoDiff<T>(scalar - dual.val_, -dual.dx_, -dual.dy_);
}

// multiplication operators
template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  return AutoDiff<T>(dual1.val_ * dual2.val_,
                     dual1.val_ * dual2.dx_ + dual1.dx_ * dual2.val_,
                     dual1.val_ * dual2.dy_ + dual1.dy_ * dual2.val_);
}

template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual, const double &scalar) {
  return AutoDiff<T>(dual.val_ * scalar, dual.dx_ * scalar, dual.dy_ * scalar);
}

template <typename T>
AutoDiff<T> operator*(const double &scalar, const AutoDiff<T> &dual) {
  return dual * scalar;
}

// division operators
template <typename T>
AutoDiff<T> operator/(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  const auto inv_val_ = 1.0 / dual2.val_;
  return AutoDiff<T>(dual1.val_ * inv_val_,
                     (dual1.dx_ - dual1.val_ * dual2.dx_ * inv_val_) * inv_val_,
                     (dual1.dy_ - dual1.val_ * dual2.dy_ * inv_val_)
                         * inv_val_);
}

template <typename T>
AutoDiff<T> operator/(const AutoDiff<T> &dual, const double &scalar) {
  const auto inv_scalar = 1.0 / scalar;
  return AutoDiff<T>(
      dual.val_ * inv_scalar, dual.dx_ * inv_scalar, dual.dy_ * inv_scalar);
}

template <typename T>
AutoDiff<T> operator/(const double &scalar, const AutoDiff<T> &dual) {
  const auto inv_val_ = 1.0 / dual.val_;
  return AutoDiff<T>(scalar * inv_val_,
                     -scalar * dual.dx_ * inv_val_ * inv_val_,
                     -scalar * dual.dy_ * inv_val_ * inv_val_);
}

// exponential function
template <typename T>
AutoDiff<T> exp(const AutoDiff<T> &x) {
  const auto exp_val_ = exp(x.val_);
  return AutoDiff<T>(exp_val_, exp_val_ * x.dx_, exp_val_ * x.dy_);
}

// square root function
template <typename T>
AutoDiff<T> sqrt(const AutoDiff<T> &x) {
  const auto sqrt_val_ = sqrt(x.val_);
  const auto inv_sqrt = 0.5 / sqrt_val_;
  return AutoDiff<T>(sqrt_val_, x.dx_ * inv_sqrt, x.dy_ * inv_sqrt);
}

// hyperbolic tangent function
template <typename T>
AutoDiff<T> tanh(const AutoDiff<T> &x) {
  const auto tanh_val_ = tanh(x.val_);
  const auto sech2_val_ = 1.0 - tanh_val_ * tanh_val_;
  return AutoDiff<T>(tanh_val_, x.dx_ * sech2_val_, x.dy_ * sech2_val_);
}

// First order derivatives
class AutoDiff1 : public AutoDiff<double> {
public:

  // Constructors
  AutoDiff1(double val_, double dx_ = 0.0, double dy_ = 0.0)
      : AutoDiff<double>(val_, dx_, dy_) {}
  AutoDiff1(const AutoDiff<double> &other)
      : AutoDiff<double>(other) {}
  // Aliases for convenient access of the results
  double &val = val_;
  double &dx = dx_;
  double &dy = dy_;
};

// Second order derivatives
class AutoDiff2 : public AutoDiff<AutoDiff1> {
public:

  // Constructors
  AutoDiff2(double val_, double dx_ = 0.0, double dy_ = 0.0)
      : AutoDiff<AutoDiff1>(
            AutoDiff1(val_, dx_, dy_), AutoDiff1(dx_), AutoDiff1(dy_)) {}
  AutoDiff2(const AutoDiff<AutoDiff1> &other)
      : AutoDiff<AutoDiff1>(other) {}
  // Aliases for convenient access of the results
  double &val = val_.val;
  double &dx = dx_.val;
  double &dy = dy_.val;
  double &dxx = dx_.dx;
  double &dxy = dx_.dy;
  double &dyy = dy_.dy;
};

#endif
