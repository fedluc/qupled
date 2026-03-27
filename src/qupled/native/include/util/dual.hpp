#ifndef DUAL_HPP
#define DUAL_HPP

#include <cmath>
#include <vector>

/**
 * @brief Recursive template class for forward-mode automatic differentiation.
 *
 * A @p Dual<Order> number carries a primal value and a gradient vector of
 * @p Dual<Order-1> entries, enabling computation of derivatives up to the
 * specified order. The base case @p Dual<1> stores a @c double primal and a
 * @c std::vector<double> gradient.
 *
 * Arithmetic operators (+, -, *, /) and common transcendental functions (exp,
 * log, sqrt, tanh) are overloaded to propagate the chain rule automatically.
 *
 * Convenience wrappers (@p Dual11, @p Dual21, @p Dual12, @p Dual22) are
 * provided for the most common combinations of order and number of variables.
 *
 * @tparam Order Differentiation order (≥ 1).
 */
template <int Order>
class Dual {
public:

  /** @brief Primal value (itself a Dual<Order-1>). */
  Dual<Order - 1> func;
  /** @brief Gradient vector (one entry per independent variable). */
  std::vector<Dual<Order - 1>> grad;

  /**
   * @brief Construct from a lower-order dual primal with a seed index.
   * @param func_  Primal value.
   * @param nvar   Number of independent variables.
   * @param index  Index of the variable seeded with 1; -1 for no seeding.
   */
  Dual(const Dual<Order - 1> &func_, const int nvar, const int index)
      : func(func_),
        grad(nvar, Dual<Order - 1>(0.0, nvar, -1)) {
    if (index >= 0 && index < nvar) {
      grad[index] = Dual<Order - 1>(1.0, nvar, -1);
    }
  }

  /**
   * @brief Construct from a scalar primal with a seed index.
   * @param func_  Scalar primal value.
   * @param nvar   Number of independent variables.
   * @param index  Index of the variable seeded with 1; -1 for no seeding.
   */
  Dual(const double &func_, const int nvar, const int index)
      : Dual(Dual<Order - 1>(func_, nvar, index), nvar, index) {}
};

/**
 * @brief Base case: first-order dual number.
 */
template <>
class Dual<1> {
public:

  /** @brief Primal (function) value. */
  double func;
  /** @brief Gradient vector over all independent variables. */
  std::vector<double> grad;

  /**
   * @brief Construct a first-order dual number.
   * @param func_  Scalar primal value.
   * @param nvar   Number of independent variables.
   * @param index  Index of the seeded variable; -1 for no seed.
   */
  Dual(const double &func_, const int nvar, const int index)
      : func(func_),
        grad(nvar) {
    if (index >= 0 && index < nvar) { grad[index] = 1; }
  }
};

// -----------------------------------------------------------------
// Arithmetic operators
// -----------------------------------------------------------------

/**
 * @brief Add two dual numbers.
 * @param dual1 Left operand.
 * @param dual2 Right operand.
 * @return Element-wise sum.
 */
template <int Order>
Dual<Order> operator+(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t nvar = dual1.grad.size();
  Dual<Order> result = Dual<Order>(dual1.func + dual2.func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual1.grad[i] + dual2.grad[i];
  }
  return result;
}

/**
 * @brief Add a dual number and a scalar.
 * @param dual   Dual operand.
 * @param scalar Scalar to add to the primal.
 * @return Dual number with @p scalar added to the primal.
 */
template <int Order>
Dual<Order> operator+(const Dual<Order> &dual, const double &scalar) {
  const size_t nvar = dual.grad.size();
  Dual<Order> result = Dual<Order>(dual.func + scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual.grad[i];
  }
  return result;
}

/**
 * @brief Add a scalar and a dual number.
 * @param scalar Scalar to add to the primal.
 * @param dual   Dual operand.
 * @return Dual number with @p scalar added to the primal.
 */
template <int Order>
Dual<Order> operator+(const double &scalar, const Dual<Order> &dual) {
  return dual + scalar;
}

/**
 * @brief Subtract @p dual2 from @p dual1.
 * @param dual1 Left operand.
 * @param dual2 Right operand.
 * @return Element-wise difference.
 */
template <int Order>
Dual<Order> operator-(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t nvar = dual1.grad.size();
  Dual<Order> result = Dual<Order>(dual1.func - dual2.func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual1.grad[i] - dual2.grad[i];
  }
  return result;
}

/**
 * @brief Subtract a scalar from a dual number.
 * @param dual   Left operand.
 * @param scalar Scalar to subtract from the primal.
 * @return Dual with @p scalar subtracted from the primal.
 */
template <int Order>
Dual<Order> operator-(const Dual<Order> &dual, double scalar) {
  const size_t nvar = dual.grad.size();
  Dual<Order> result = Dual<Order>(dual.func - scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual.grad[i];
  }
  return result;
}

/**
 * @brief Subtract a dual number from a scalar.
 * @param scalar Left operand.
 * @param dual   Right operand.
 * @return Negated dual with @p scalar added to the primal.
 */
template <int Order>
Dual<Order> operator-(const double &scalar, const Dual<Order> &dual) {
  const size_t nvar = dual.grad.size();
  Dual<Order> result = Dual<Order>(scalar - dual.func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = -1.0 * dual.grad[i];
  }
  return result;
}

/**
 * @brief Multiply two dual numbers using the product rule.
 * @param dual1 Left operand.
 * @param dual2 Right operand.
 * @return Product dual number.
 */
template <int Order>
Dual<Order> operator*(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t nvar = dual1.grad.size();
  Dual<Order> result = Dual<Order>(dual1.func * dual2.func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual1.func * dual2.grad[i] + dual1.grad[i] * dual2.func;
  }
  return result;
}

/**
 * @brief Multiply a dual number by a scalar.
 * @param dual   Dual operand.
 * @param scalar Scalar multiplier.
 * @return Scaled dual number.
 */
template <int Order>
Dual<Order> operator*(const Dual<Order> &dual, const double &scalar) {
  const size_t nvar = dual.grad.size();
  Dual<Order> result = Dual<Order>(dual.func * scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual.grad[i] * scalar;
  }
  return result;
}

/**
 * @brief Multiply a scalar by a dual number.
 * @param scalar Scalar multiplier.
 * @param dual   Dual operand.
 * @return Scaled dual number.
 */
template <int Order>
Dual<Order> operator*(const double &scalar, const Dual<Order> &dual) {
  return dual * scalar;
}

/**
 * @brief Divide @p dual1 by @p dual2 using the quotient rule.
 * @param dual1 Numerator.
 * @param dual2 Denominator.
 * @return Quotient dual number.
 */
template <int Order>
Dual<Order> operator/(const Dual<Order> &dual1, const Dual<Order> &dual2) {
  const size_t nvar = dual1.grad.size();
  const auto inv_func = 1.0 / dual2.func;
  Dual<Order> result = Dual<Order>(dual1.func * inv_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] =
        (dual1.grad[i] - dual1.func * dual2.grad[i] * inv_func) * inv_func;
  }
  return result;
}

/**
 * @brief Divide a dual number by a scalar.
 * @param dual   Numerator.
 * @param scalar Denominator scalar.
 * @return Scaled dual number.
 */
template <int Order>
Dual<Order> operator/(const Dual<Order> &dual, const double &scalar) {
  const size_t nvar = dual.grad.size();
  const auto inv_scalar = 1.0 / scalar;
  Dual<Order> result = Dual<Order>(dual.func * inv_scalar, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = dual.grad[i] * inv_scalar;
  }
  return result;
}

/**
 * @brief Divide a scalar by a dual number.
 * @param scalar Numerator scalar.
 * @param dual   Denominator dual number.
 * @return Quotient dual number.
 */
template <int Order>
Dual<Order> operator/(const double &scalar, const Dual<Order> &dual) {
  const size_t nvar = dual.grad.size();
  const auto inv_func = 1.0 / dual.func;
  const auto inv_func2 = inv_func * inv_func;
  Dual<Order> result = Dual<Order>(scalar * inv_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = -scalar * dual.grad[i] * inv_func2;
  }
  return result;
}

// -----------------------------------------------------------------
// Transcendental functions
// -----------------------------------------------------------------

/**
 * @brief Dual exponential @f$e^x@f$.
 * @param x Argument.
 * @return Dual number with primal @f$e^{x.func}@f$ and propagated gradient.
 */
template <int Order>
Dual<Order> exp(const Dual<Order> &x) {
  const size_t nvar = x.grad.size();
  const auto exp_func = exp(x.func);
  Dual<Order> result = Dual<Order>(exp_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = exp_func * x.grad[i];
  }
  return result;
}

/**
 * @brief Dual natural logarithm @f$\ln x@f$.
 * @param x Argument.
 * @return Dual number with primal @f$\ln(x.func)@f$ and propagated gradient.
 */
template <int Order>
Dual<Order> log(const Dual<Order> &x) {
  const size_t nvar = x.grad.size();
  const auto log_func = log(x.func);
  Dual<Order> result = Dual<Order>(log_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = x.grad[i] / x.func;
  }
  return result;
}

/**
 * @brief Dual square root @f$\sqrt{x}@f$.
 * @param x Argument.
 * @return Dual number with primal @f$\sqrt{x.func}@f$ and propagated gradient.
 */
template <int Order>
Dual<Order> sqrt(const Dual<Order> &x) {
  const size_t nvar = x.grad.size();
  const auto sqrt_func = sqrt(x.func);
  const auto inv_sqrt = 0.5 / sqrt_func;
  Dual<Order> result = Dual<Order>(sqrt_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = x.grad[i] * inv_sqrt;
  }
  return result;
}

/**
 * @brief Dual hyperbolic tangent @f$\tanh x@f$.
 * @param x Argument.
 * @return Dual number with primal @f$\tanh(x.func)@f$ and propagated gradient.
 */
template <int Order>
Dual<Order> tanh(const Dual<Order> &x) {
  const size_t nvar = x.grad.size();
  const auto tanh_func = tanh(x.func);
  const auto sech2_func = 1.0 - tanh_func * tanh_func;
  Dual<Order> result = Dual<Order>(tanh_func, nvar, -1);
  for (size_t i = 0; i < nvar; ++i) {
    result.grad[i] = x.grad[i] * sech2_func;
  }
  return result;
}

// -----------------------------------------------------------------
// Convenience wrappers for common (order, nvar) combinations
// -----------------------------------------------------------------

/**
 * @brief First-order dual number for a function of one variable.
 *
 * Provides named accessors @p val() and @p dx().
 */
class Dual11 : public Dual<1> {
public:

  /**
   * @brief Construct a Dual11 number.
   * @param func_  Scalar primal value.
   * @param index  0 to seed @f$dx@f$; -1 for a constant.
   */
  explicit Dual11(const double func_, const int index = -1)
      : Dual<1>(func_, 1, index) {}

  /**
   * @brief Construct from a generic @p Dual<1>.
   * @param other Source dual number.
   */
  Dual11(const Dual<1> &other)
      : Dual<1>(other) {}

  /** @brief Return the primal value @f$f@f$. */
  const double &val() const { return func; }
  /** @brief Return the first derivative @f$\partial f / \partial x@f$. */
  const double &dx() const { return grad[0]; }
};

/**
 * @brief Second-order dual number for a function of one variable.
 *
 * Provides named accessors @p val(), @p dx(), and @p dxx().
 */
class Dual21 : public Dual<2> {
public:

  /**
   * @brief Construct a Dual21 number.
   * @param func_  Scalar primal value.
   * @param index  0 to seed @f$dx@f$; -1 for a constant.
   */
  explicit Dual21(const double func_, const int index = -1)
      : Dual<2>(func_, 1, index) {}

  /**
   * @brief Construct from a generic @p Dual<2>.
   * @param other Source dual number.
   */
  Dual21(const Dual<2> &other)
      : Dual<2>(other) {}

  /** @brief Return the primal value @f$f@f$. */
  const double &val() const { return func.func; }
  /** @brief Return the first derivative @f$\partial f / \partial x@f$. */
  const double &dx() const { return grad[0].func; }
  /** @brief Return the second derivative @f$\partial^2 f / \partial x^2@f$. */
  const double &dxx() const { return grad[0].grad[0]; }
};

/**
 * @brief First-order dual number for a function of two variables.
 *
 * Provides named accessors @p val(), @p dx(), and @p dy().
 */
class Dual12 : public Dual<1> {
public:

  /**
   * @brief Construct a Dual12 number.
   * @param func_  Scalar primal value.
   * @param index  0 to seed @f$dx@f$, 1 to seed @f$dy@f$; -1 for a constant.
   */
  explicit Dual12(const double func_, const int index = -1)
      : Dual<1>(func_, 2, index) {}

  /**
   * @brief Construct from a generic @p Dual<1>.
   * @param other Source dual number.
   */
  Dual12(const Dual<1> &other)
      : Dual<1>(other) {}

  /** @brief Return the primal value @f$f@f$. */
  const double &val() const { return func; }
  /** @brief Return @f$\partial f / \partial x@f$. */
  const double &dx() const { return grad[0]; }
  /** @brief Return @f$\partial f / \partial y@f$. */
  const double &dy() const { return grad[1]; }
};

/**
 * @brief Second-order dual number for a function of two variables.
 *
 * Provides named accessors @p val(), @p dx(), @p dy(), @p dxx(), @p dxy(),
 * @p dyx(), and @p dyy().
 */
class Dual22 : public Dual<2> {
public:

  /**
   * @brief Construct a Dual22 number.
   * @param func_  Scalar primal value.
   * @param index  0 to seed @f$dx@f$, 1 to seed @f$dy@f$; -1 for a constant.
   */
  explicit Dual22(const double &func_, const int index = -1)
      : Dual<2>(func_, 2, index) {}

  /**
   * @brief Construct from a generic @p Dual<2>.
   * @param other Source dual number.
   */
  Dual22(const Dual<2> &other)
      : Dual<2>(other) {}

  /** @brief Return the primal value @f$f@f$. */
  const double &val() const { return func.func; }
  /** @brief Return @f$\partial f / \partial x@f$. */
  const double &dx() const { return grad[0].func; }
  /** @brief Return @f$\partial f / \partial y@f$. */
  const double &dy() const { return grad[1].func; }
  /** @brief Return @f$\partial^2 f / \partial x^2@f$. */
  const double &dxx() const { return grad[0].grad[0]; }
  /** @brief Return @f$\partial^2 f / \partial x \partial y@f$. */
  const double &dxy() const { return grad[0].grad[1]; }
  /** @brief Return @f$\partial^2 f / \partial y \partial x@f$. */
  const double &dyx() const { return grad[1].grad[0]; }
  /** @brief Return @f$\partial^2 f / \partial y^2@f$. */
  const double &dyy() const { return grad[1].grad[1]; }
};

#endif
