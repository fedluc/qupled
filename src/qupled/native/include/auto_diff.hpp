#ifndef AUTO_DIFF_HPP
#define AUTO_DIFF_HPP

#include <cmath>
#include <vector>

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

template <int Order>
class Dual {
public:

  double val;
  std::vector<Dual<Order-1>> der;
  // Constructors
  Dual(const double &val_, const size_t &nvar, const size_t &index)
      : val(val_),
        der(nvar, Dual<Order-1>(0.0, nvar, -1)) {
    if (index >= 0 && index < nvar) {
            der[index] = Dual<Order - 1>(1.0, nvar, -1);
    }
  }
};

template <>
class Dual<1> {
public:

  double val;
  std::vector<double> der;
  // Constructors
  Dual(const double &val_, const size_t &nvar, const size_t &index)
      : val(val_),
        der(nvar) {
    if (index >= 0 && index < nvar) { der[index] = 1; }
  }
};


// addition operators
template <int Order>
Dual<Order> operator+(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t &nvar = dual1.der.size();
  Dual<Order> result = Dual<Order>(dual1.val + dual2.val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual1.der[i] + dual2.der[i];
  }
  return result;
}

template <int Order>
Dual<Order> operator+(const Dual<Order> &dual, const double &scalar) {
  const size_t &nvar = dual.der.size();
  Dual<Order> result = Dual<Order>(dual.val + scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual.der[i];
  }
  return result;
}

template <int Order>
Dual<Order> operator+(const double &scalar, const Dual<Order> &dual) {
  return dual + scalar;
}

// subtraction operators
template <int Order>
Dual<Order> operator-(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t &nvar = dual1.der.size();
  Dual<Order> result = Dual<Order>(dual1.val - dual2.val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual1.der[i] - dual2.der[i];
  }
  return result;
}

template <int Order>
Dual<Order> operator-(const Dual<Order> &dual, double scalar) {
  const size_t &nvar = dual.der.size();
  Dual<Order> result = Dual<Order>(dual.val - scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual.der[i];
  }
  return result;
}

template <int Order>
Dual<Order> operator-(const double &scalar, const Dual<Order> &dual) {
  const size_t &nvar = dual.der.size();
  Dual<Order> result = Dual<Order>(scalar - dual.val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = -dual.der[i];
  }
  return result;
}

// multiplication operators
template <int Order>
Dual<Order> operator*(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t &nvar = dual1.der.size();
  Dual<Order> result = Dual<Order>(dual1.val * dual2.val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual1.val * dual2.der[i] + dual1.der[i] * dual2.val;
  }
  return result;
}

template <int Order>
Dual<Order> operator*(const Dual<Order> &dual, const double &scalar) {
  const size_t &nvar = dual.der.size();
  Dual<Order> result = Dual<Order>(dual.val * scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual.der[i] * scalar;
  }
  return result;
}

template <int Order>
Dual<Order> operator*(const double &scalar, const Dual<Order> &dual) {
  return dual * scalar;
}

// division operators
template <int Order>
Dual<Order> operator/(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t &nvar = dual1.der.size();
  const double inv_val = 1.0 / dual2.val;
  Dual<Order> result = Dual<Order>(dual1.val * inv_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] =
        (dual1.der[i] - dual1.val * dual2.der[i] * inv_val) * inv_val;
  }
  return result;
}

template <int Order>
Dual<Order> operator/(const Dual<Order> &dual, const double &scalar) {
  const size_t &nvar = dual.der.size();
  const double inv_scalar = 1.0 / scalar;
  Dual<Order> result = Dual<Order>(dual.val * inv_scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = dual.der[i] * inv_scalar;
  }
  return result;
}

template <int Order>
Dual<Order> operator/(const double &scalar, const Dual<Order> &dual) {
  const size_t &nvar = dual.der.size();
  const double inv_val = 1.0 / dual.val;
  const double inv_val2 = inv_val * inv_val;
  Dual<Order> result = Dual<Order>(scalar * inv_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = -scalar * dual.der[i] * inv_val2;
  }
  return result;
}

// exponential function
template <int Order>
Dual<Order> exp(const Dual<Order> &x) {
  const size_t &nvar = x.der.size();
  const double exp_val = exp(x.val);
  Dual<Order> result = Dual<Order>(exp_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = exp_val * x.der[i];
  }
  return result;
}

// square root function
template <int Order>
Dual<Order> sqrt(const Dual<Order> &x) {
  const size_t &nvar = x.der.size();
  const double sqrt_val = sqrt(x.val);
  const double inv_sqrt = 0.5 / sqrt_val;
  Dual<Order> result = Dual<Order>(sqrt_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = x.der[i] * inv_sqrt;
  }
  return result;
}

// hyperbolic tangent function
template <int Order>
Dual<Order> tanh(const Dual<Order> &x) {
  const size_t &nvar = x.der.size();
  const double tanh_val = tanh(x.val);
  const double sech2_val = 1.0 - tanh_val * tanh_val;
  Dual<Order> result = Dual<Order>(tanh_val, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.der[i] = x.der[i] * sech2_val;
  }
  return result;
}

// First order derivatives
class AutoDiff1 : public Dual<1> {
public:

  // Constructors
  AutoDiff1(const double val_, const size_t index = -1)
      : Dual<1>(val_, 2, index) {}
  AutoDiff1(const Dual<1> &other)
      : Dual<1>(other) {}
  // Aliases for convenient access of the results
  const double &dx() const { return der[0]; }
  const double &dy() const { return der[1]; }
};

// Second order derivatives
class AutoDiff2 : public Dual<2> {
public:

  // Constructors
  AutoDiff2(const double &val_, const size_t index = -1)
      : Dual<2>(val_, 2, index) {}
  AutoDiff2(const Dual<2> &other)
      : Dual<2>(other) {}
  // Aliases for convenient access of the results
  const double &dx() const { return der[0].val; }
  const double &dy() const { return der[1].val; }
  const double &dxx() const { return der[0].der[0]; }
  const double &dxy() const { return der[0].der[1]; }
  const double &dyx() const { return der[1].der[0]; }
  const double &dyy() const { return der[1].der[1]; }
};

#endif
