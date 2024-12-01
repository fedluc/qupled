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

  double val;
  std::vector<T> der;
  // Constructors
  AutoDiff(const double &val_, const size_t &nvar, const size_t &index)
      : val(val_),
        der(nvar) {
    if (index >= 0 && index < nvar) { der[index] = 1; }
  }
};

// addition operators
template <typename T>
AutoDiff<T> operator+(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  const size_t &nvar = dual1.der.size();
  AutoDiff<T> result = AutoDiff<T>(dual1.val + dual2.val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual1.der[i] + dual2.der[i];
  }
  return result;
}

template <typename T>
AutoDiff<T> operator+(const AutoDiff<T> &dual, const double &scalar) {
  const size_t &nvar = dual.der.size();
  AutoDiff<T> result = AutoDiff<T>(dual.val + scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
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
  const size_t &nvar = dual1.der.size();
  AutoDiff<T> result = AutoDiff<T>(dual1.val - dual2.val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual1.der[i] - dual2.der[i];
  }
  return result;
}

template <typename T>
AutoDiff<T> operator-(const AutoDiff<T> &dual, double scalar) {
  const size_t &nvar = dual.der.size();
  AutoDiff<T> result = AutoDiff<T>(dual.val - scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual.der[i];
  }
  return result;
}

template <typename T>
AutoDiff<T> operator-(const double &scalar, const AutoDiff<T> &dual) {
  const size_t &nvar = dual.der.size();
  AutoDiff<T> result = AutoDiff<T>(scalar - dual.val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = -dual.der[i];
  }
  return result;
}

// multiplication operators
template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual1, const AutoDiff<T> &dual2) {
  const size_t &nvar = dual1.der.size();
  AutoDiff<T> result = AutoDiff<T>(dual1.val * dual2.val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual1.val * dual2.der[i] + dual1.der[i] * dual2.val;
  }
  return result;
}

template <typename T>
AutoDiff<T> operator*(const AutoDiff<T> &dual, const double &scalar) {
  const size_t &nvar = dual.der.size();
  AutoDiff<T> result = AutoDiff<T>(dual.val * scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
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
  const size_t &nvar = dual1.der.size();
  const double inv_val = 1.0 / dual2.val;
  AutoDiff<T> result = AutoDiff<T>(dual1.val * inv_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] =
        (dual1.der[i] - dual1.val * dual2.der[i] * inv_val) * inv_val;
  }
  return result;
}

template <typename T>
AutoDiff<T> operator/(const AutoDiff<T> &dual, const double &scalar) {
  const size_t &nvar = dual.der.size();
  const double inv_scalar = 1.0 / scalar;
  AutoDiff<T> result = AutoDiff<T>(dual.val * inv_scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual.der[i] * inv_scalar;
  }
  return result;
}

template <typename T>
AutoDiff<T> operator/(const double &scalar, const AutoDiff<T> &dual) {
  const size_t &nvar = dual.der.size();
  const double inv_val = 1.0 / dual.val;
  const double inv_val2 = inv_val * inv_val;
  AutoDiff<T> result = AutoDiff<T>(scalar * inv_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = -scalar * dual.der[i] * inv_val2;
  }
  return result;
}

// exponential function
template <typename T>
AutoDiff<T> exp(const AutoDiff<T> &x) {
  const size_t &nvar = x.der.size();
  const double exp_val = exp(x.val);
  AutoDiff<T> result = AutoDiff<T>(exp_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = exp_val * x.der[i];
  }
  return result;
}

// square root function
template <typename T>
AutoDiff<T> sqrt(const AutoDiff<T> &x) {
  const size_t &nvar = x.der.size();
  const double sqrt_val = sqrt(x.val);
  const double inv_sqrt = 0.5 / sqrt_val;
  AutoDiff<T> result = AutoDiff<T>(sqrt_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = x.der[i] * inv_sqrt;
  }
  return result;
}

// hyperbolic tangent function
template <typename T>
AutoDiff<T> tanh(const AutoDiff<T> &x) {
  const size_t &nvar = x.der.size();
  const double tanh_val = tanh(x.val);
  const double sech2_val = 1.0 - tanh_val * tanh_val;
  AutoDiff<T> result = AutoDiff<T>(tanh_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = x.der[i] * sech2_val;
  }
  return result;
}

// First order derivatives
class AutoDiff1 : public AutoDiff<double> {
public:

  // Constructors
  AutoDiff1(const double val_ = 0, const size_t index = -1)
      : AutoDiff<double>(val_, 2, index) {}
  AutoDiff1(const AutoDiff<double> &other)
      : AutoDiff<double>(other) {}
  // Aliases for convenient access of the results
  const double &dx() const { return der[0]; }
  const double &dy() const { return der[1]; }
};

// Second order derivatives
class AutoDiff2 : public AutoDiff<AutoDiff1> {
public:

  // Constructors
  AutoDiff2(const double &val_, const size_t index = -1)
      : AutoDiff<AutoDiff1>(val_, 2, index) {}
  AutoDiff2(const AutoDiff<AutoDiff1> &other)
      : AutoDiff<AutoDiff1>(other) {}
  // Aliases for convenient access of the results
  const double &dx() const { return der[0].val; }
  const double &dy() const { return der[1].val; }
  const double &dxx() const { return der[0].der[0]; }
  const double &dxy() const { return der[0].der[1]; }
  const double &dyx() const { return der[1].der[0]; }
  const double &dyy() const { return der[1].der[1]; }
};

#endif
