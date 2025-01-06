#ifndef DUAL_HPP
#define DUAL_HPP

#include "numerics.hpp"
#include <cmath>
#include <vector>

using namespace SpecialFunctions;

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

template <int Order>
class Dual {
public:

  Dual<Order - 1> func;
  std::vector<Dual<Order - 1>> grad;
  // Constructors
  Dual(const Dual<Order - 1> &func_, const int nvar, const int index)
      : func(func_),
        grad(nvar, Dual<Order - 1>(0.0, nvar, -1)) {
    if (index >= 0 && index < nvar) {
      grad[index] = Dual<Order - 1>(1.0, nvar, -1);
    }
  }
  Dual(const double &func_, const int nvar, const int index)
      : Dual(Dual<Order - 1>(func_, nvar, index), nvar, index) {}
};

template <>
class Dual<1> {
public:

  double func;
  std::vector<double> grad;
  // Constructors
  Dual(const double &func_, const int nvar, const int index)
      : func(func_),
        grad(nvar) {
    if (index >= 0 && index < nvar) { grad[index] = 1; }
  }
};

// addition operators
template <int Order>
Dual<Order> operator+(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t &nvar = dual1.grad.size();
  Dual<Order> result = Dual<Order>(dual1.func + dual2.func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual1.grad[i] + dual2.grad[i];
  }
  return result;
}

template <int Order>
Dual<Order> operator+(const Dual<Order> &dual, const double &scalar) {
  const size_t &nvar = dual.grad.size();
  Dual<Order> result = Dual<Order>(dual.func + scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual.grad[i];
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
  const size_t &nvar = dual1.grad.size();
  Dual<Order> result = Dual<Order>(dual1.func - dual2.func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual1.grad[i] - dual2.grad[i];
  }
  return result;
}

template <int Order>
Dual<Order> operator-(const Dual<Order> &dual, double scalar) {
  const size_t &nvar = dual.grad.size();
  Dual<Order> result = Dual<Order>(dual.func - scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual.grad[i];
  }
  return result;
}

template <int Order>
Dual<Order> operator-(const double &scalar, const Dual<Order> &dual) {
  const size_t &nvar = dual.grad.size();
  Dual<Order> result = Dual<Order>(scalar - dual.func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = -1.0 * dual.grad[i];
  }
  return result;
}

// multiplication operators
template <int Order>
Dual<Order> operator*(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t &nvar = dual1.grad.size();
  Dual<Order> result = Dual<Order>(dual1.func * dual2.func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual1.func * dual2.grad[i] + dual1.grad[i] * dual2.func;
  }
  return result;
}

template <int Order>
Dual<Order> operator*(const Dual<Order> &dual, const double &scalar) {
  const size_t &nvar = dual.grad.size();
  Dual<Order> result = Dual<Order>(dual.func * scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual.grad[i] * scalar;
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
  const size_t &nvar = dual1.grad.size();
  const auto &inv_func = 1.0 / dual2.func;
  Dual<Order> result = Dual<Order>(dual1.func * inv_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] =
        (dual1.grad[i] - dual1.func * dual2.grad[i] * inv_func) * inv_func;
  }
  return result;
}

template <int Order>
Dual<Order> operator/(const Dual<Order> &dual, const double &scalar) {
  const size_t &nvar = dual.grad.size();
  const auto &inv_scalar = 1.0 / scalar;
  Dual<Order> result = Dual<Order>(dual.func * inv_scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual.grad[i] * inv_scalar;
  }
  return result;
}

template <int Order>
Dual<Order> operator/(const double &scalar, const Dual<Order> &dual) {
  const size_t &nvar = dual.grad.size();
  const auto &inv_func = 1.0 / dual.func;
  const auto &inv_func2 = inv_func * inv_func;
  Dual<Order> result = Dual<Order>(scalar * inv_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = -scalar * dual.grad[i] * inv_func2;
  }
  return result;
}

// Smaller-than operator
template <int Order>
std::vector<bool> operator<(const Dual<Order> &dual1,
                            const Dual<Order> &dual2) {
  const size_t &nvar = dual1.grad.size();
  std::vector<bool> result(nvar + 1);
  result[0] = dual1.func < dual2.func;
  for (size_t i = 1; i < nvar + 1; ++i) {
    result[i] = dual1.grad[i] < dual1.grad[i];
  }
  return result;
}

template <int Order>
std::vector<bool> operator<(const Dual<Order> &dual, const double &scalar) {
  const auto scalarDual = Dual<Order>(scalar, dual.grad.size(), -1);
  return dual < scalarDual;
}

template <int Order>
std::vector<bool> operator<(const double &scalar, const Dual<Order> &dual) {
  return dual > scalar;
}

// Larger-than operator
template <int Order>
std::vector<bool> operator>(const Dual<Order> &dual1,
                            const Dual<Order> &dual2) {
  const size_t &nvar = dual1.grad.size();
  std::vector<bool> result(nvar + 1);
  result[0] = dual1.func > dual2.func;
  for (size_t i = 1; i < nvar + 1; ++i) {
    result[i] = dual1.grad[i] > dual1.grad[i];
  }
  return result;
}

template <int Order>
std::vector<bool> operator>(const Dual<Order> &dual, const double &scalar) {
  const auto scalarDual = Dual<Order>(scalar, dual.grad.size(), -1);
  return dual > scalarDual;
}

template <int Order>
std::vector<bool> operator>(const double &scalar, const Dual<Order> &dual) {
  return dual < scalar;
}

// exponential function
template <int Order>
Dual<Order> exp(const Dual<Order> &x) {
  const size_t &nvar = x.grad.size();
  const auto &exp_func = exp(x.func);
  Dual<Order> result = Dual<Order>(exp_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = exp_func * x.grad[i];
  }
  return result;
}

// logarithmic function
template <int Order>
Dual<Order> log(const Dual<Order> &x) {
  const size_t &nvar = x.grad.size();
  const auto &log_func = log(x.func);
  Dual<Order> result = Dual<Order>(log_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = x.grad[i] / x.func;
  }
  return result;
}

// square root function
template <int Order>
Dual<Order> sqrt(const Dual<Order> &x) {
  const size_t &nvar = x.grad.size();
  const auto &sqrt_func = sqrt(x.func);
  const auto &inv_sqrt = 0.5 / sqrt_func;
  Dual<Order> result = Dual<Order>(sqrt_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = x.grad[i] * inv_sqrt;
  }
  return result;
}

// hyperbolic tangent function
template <int Order>
Dual<Order> tanh(const Dual<Order> &x) {
  const size_t &nvar = x.grad.size();
  const auto &tanh_func = tanh(x.func);
  const auto &sech2_func = 1.0 - tanh_func * tanh_func;
  Dual<Order> result = Dual<Order>(tanh_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = x.grad[i] * sech2_func;
  }
  return result;
}

// Absolute value function
template <int Order>
Dual<Order> abs(const Dual<Order> &x) {
  const size_t &nvar = x.grad.size();
  const auto abs_func = abs(x.func);
  const auto sign_func = x.func / abs_func;
  Dual<Order> result = Dual<Order>(abs_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = x.grad[i] * sign_func;
  }
  return result;
}

// dilog function
template <int Order>
Dual<Order> dilog(const Dual<Order> &x) {
  const size_t nvar = x.grad.size();
  const auto func_value = dilog(x.func);
  const auto derivative_value = -1.0 * log(1.0 - x.func) / x.func;
  Dual<Order> result(func_value, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = x.grad[i] * derivative_value;
  }
  return result;
}

// spence function
template <int Order>
Dual<Order> spence(const Dual<Order> &x) {
  const bool smallerThanOne = (x.func < 1.0)[0];
  if (smallerThanOne) {
    Dual<Order> result = dilog(x);
    return result;
  }
  const auto logx = log(x);
  const auto logx2 = logx * logx;
  const auto inv_x = 1.0 / x;
  const auto dilog_inv_x = dilog(inv_x);
  return M_PI * M_PI / 3.0 - 0.5 * logx2 - dilog_inv_x;
}

// -----------------------------------------------------------------
// Wrappers for specific derivative calculations
// -----------------------------------------------------------------

// First order derivatives for functions of one variable
class Dual11 : public Dual<1> {
public:

  // Constructors
  explicit Dual11(const double func_, const int index = -1)
      : Dual<1>(func_, 1, index) {}
  Dual11(const Dual<1> &other)
      : Dual<1>(other) {}
  // Aliases for convenient access of the results
  const double &val() const { return func; }
  const double &dx() const { return grad[0]; }
};

// Second order derivatives for functions of one variable
class Dual21 : public Dual<2> {
public:

  // Constructors
  explicit Dual21(const double func_, const int index = -1)
      : Dual<2>(func_, 1, index) {}
  Dual21(const Dual<2> &other)
      : Dual<2>(other) {}
  Dual21(const double val, const double dx, const double dxx)
      : Dual21(0) {
    func.func = val;
    grad[0].func = dx;
    grad[0].grad[0] = dxx;
  }
  // Aliases for convenient access of the results
  const double &val() const { return func.func; }
  const double &dx() const { return grad[0].func; }
  const double &dxx() const { return grad[0].grad[0]; }
};

// First order derivatives for functions of two variables
class Dual12 : public Dual<1> {
public:

  // Constructors
  explicit Dual12(const double func_, const int index = -1)
      : Dual<1>(func_, 2, index) {}
  Dual12(const Dual<1> &other)
      : Dual<1>(other) {}
  // Aliases for convenient access of the results
  const double &val() const { return func; }
  const double &dx() const { return grad[0]; }
  const double &dy() const { return grad[1]; }
};

// Second order derivatives for functions of two variables
class Dual22 : public Dual<2> {
public:

  // Constructors
  explicit Dual22(const double &func_, const int index = -1)
      : Dual<2>(func_, 2, index) {}
  Dual22(const Dual<2> &other)
      : Dual<2>(other) {}
  // Aliases for convenient access of the results
  const double &val() const { return func.func; }
  const double &dx() const { return grad[0].func; }
  const double &dy() const { return grad[1].func; }
  const double &dxx() const { return grad[0].grad[0]; }
  const double &dxy() const { return grad[0].grad[1]; }
  const double &dyx() const { return grad[1].grad[0]; }
  const double &dyy() const { return grad[1].grad[1]; }
};

#endif
