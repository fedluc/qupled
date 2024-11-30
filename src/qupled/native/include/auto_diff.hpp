#ifndef AUTO_DIFF_HPP
#define AUTO_DIFF_HPP

#include <vector>
#include <cmath>

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

// Automatic differentiation for two dimensional functions
class AutoDiff2 {
public:

  double val;
  double dx;
  double dy;
  double dxx;
  double dyy;
  double dxy;

  // Constructors
  AutoDiff2(double val_,
            double dx_ = 0.0,
            double dy_ = 0.0,
            double dxx_ = 0.0,
            double dyy_ = 0.0,
            double dxy_ = 0.0)
      : val(val_),
        dx(dx_),
        dy(dy_),
        dxx(dxx_),
        dyy(dyy_),
        dxy(dxy_) {}
};

// addition operators
AutoDiff2 operator+(const AutoDiff2 &dual1, const AutoDiff2 &dual2);
AutoDiff2 operator+(const AutoDiff2 &dual, const double &scalar);
AutoDiff2 operator+(const double &scalar, const AutoDiff2 &dual);

// subtraction operators
AutoDiff2 operator-(const AutoDiff2 &dual1, const AutoDiff2 &dual2);
AutoDiff2 operator-(const AutoDiff2 &dual, double scalar);
AutoDiff2 operator-(const double &scalar, const AutoDiff2 &dual);

// multiplication operators
AutoDiff2 operator*(const AutoDiff2 &dual1, const AutoDiff2 &dual2);
AutoDiff2 operator*(const AutoDiff2 &dual, const double &scalar);
AutoDiff2 operator*(const double &scalar, const AutoDiff2 &dual);

// division operators
AutoDiff2 operator/(const AutoDiff2 &dual1, const AutoDiff2 &dual2);
AutoDiff2 operator/(const AutoDiff2 &dual, const double &scalar);
AutoDiff2 operator/(const double &scalar, const AutoDiff2 &dual);

// exponential function
AutoDiff2 exp(const AutoDiff2 &x);

// square root function
AutoDiff2 sqrt(const AutoDiff2 &x);

// hyperbolic tangent function
AutoDiff2 tanh(const AutoDiff2 &x);


// -----------------------------------------------------------------
// Classes for automatic differentiation (MORE GENERAL SOLUTION
// -----------------------------------------------------------------

// Recursive Dual structure for higher-order der
template <int Order>
struct Dual {
    double val;
  std::vector<Dual<Order - 1>> der;
    // Constructor
    Dual(double v, int n, int index) : val(v), der(n, Dual<Order - 1>(0.0, n, -1)) {
        if (index >= 0 && index < n) {
            der[index] = Dual<Order - 1>(1.0, n, index); // Set derivative to 1 for the corresponding variable
        }
    }
    // Base Constructor for scalar vals
    Dual(double v) : val(v), der() {}
};

// Recursive Dual structure for higher-order der
template <>
struct Dual<1> {
  double val;
  std::vector<double> der;
    // Constructor
    Dual(double v, int n, int index) : val(v), der(n, 0.0) {
        if (index >= 0 && index < n) {
	  der[index] = 1.0; // Set derivative to 1 for the corresponding variable
        }
    }
    // Base Constructor for scalar vals
    Dual(double v) : val(v), der() {}
};

// Addition operators
template <int Order>
Dual<Order> operator+(const Dual<Order>& dual1, const Dual<Order>& dual2) {
    Dual<Order> result(dual1.val + dual2.val, dual1.der.size(), -1);
    for (size_t i = 0; i < dual1.der.size(); ++i) {
        result.der[i] = dual1.der[i] + dual2.der[i];
    }
    return result;
}

template <int Order>
Dual<Order> operator+(const Dual<Order>& dual, const double& scalar) {
    Dual<Order> result(dual.val + scalar, dual.der.size(), -1);
    result.der = dual.der;
    return result;
}

template <int Order>
Dual<Order> operator+(const double& scalar, const Dual<Order>& dual) {
    return dual + scalar;
}

// Subtraction operators
template <int Order>
Dual<Order> operator-(const Dual<Order>& dual1, const Dual<Order>& dual2) {
    Dual<Order> result(dual1.val - dual2.val, dual1.der.size(), -1);
    for (size_t i = 0; i < dual1.der.size(); ++i) {
        result.der[i] = dual1.der[i] - dual2.der[i];
    }
    return result;
}

template <int Order>
Dual<Order> operator-(const Dual<Order>& dual, const double& scalar) {
    Dual<Order> result(dual.val - scalar, dual.der.size(), -1);
    result.der = dual.der;
    return result;
}

template <int Order>
Dual<Order> operator-(const double& scalar, const Dual<Order>& dual) {
    Dual<Order> result(scalar - dual.val, dual.der.size(), -1);
    for (size_t i = 0; i < dual.der.size(); ++i) {
        result.der[i] = -dual.der[i];
    }
    return result;
}

// Multiplication operators
template <int Order>
Dual<Order> operator*(const Dual<Order>& dual1, const Dual<Order>& dual2) {
    Dual<Order> result(dual1.val * dual2.val, dual1.der.size(), -1);
    for (size_t i = 0; i < dual1.der.size(); ++i) {
        result.der[i] = dual1.der[i] * dual2.val + dual1.val * dual2.der[i];
    }
    return result;
}

template <int Order>
Dual<Order> operator*(const Dual<Order>& dual, const double& scalar) {
    Dual<Order> result(dual.val * scalar, dual.der.size(), -1);
    for (size_t i = 0; i < dual.der.size(); ++i) {
        result.der[i] = dual.der[i] * scalar;
    }
    return result;
}

template <int Order>
Dual<Order> operator*(const double& scalar, const Dual<Order>& dual) {
    return dual * scalar;
}

// Division operators
template <int Order>
Dual<Order> operator/(const Dual<Order>& dual1, const Dual<Order>& dual2) {
    Dual<Order> result(dual1.val / dual2.val, dual1.der.size(), -1);
    for (size_t i = 0; i < dual1.der.size(); ++i) {
        result.der[i] = (dual1.der[i] * dual2.val - dual1.val * dual2.der[i]) / (dual2.val * dual2.val);
    }
    return result;
}

template <int Order>
Dual<Order> operator/(const Dual<Order>& dual, const double& scalar) {
    Dual<Order> result(dual.val / scalar, dual.der.size(), -1);
    for (size_t i = 0; i < dual.der.size(); ++i) {
        result.der[i] = dual.der[i] / scalar;
    }
    return result;
}

template <int Order>
Dual<Order> operator/(const double& scalar, const Dual<Order>& dual) {
    Dual<Order> result(scalar / dual.val, dual.der.size(), -1);
    for (size_t i = 0; i < dual.der.size(); ++i) {
        result.der[i] = -scalar * dual.der[i] / (dual.val * dual.val);
    }
    return result;
}

// Exponential function
template <int Order>
Dual<Order> exp(const Dual<Order>& x) {
    Dual<Order> result(std::exp(x.val), x.der.size(), -1);
    for (size_t i = 0; i < x.der.size(); ++i) {
        result.der[i] = result.val * x.der[i];
    }
    return result;
}

// Square root function
template <int Order>
Dual<Order> sqrt(const Dual<Order>& x) {
    Dual<Order> result(std::sqrt(x.val), x.der.size(), -1);
    for (size_t i = 0; i < x.der.size(); ++i) {
        result.der[i] = x.der[i] / (2 * result.val);
    }
    return result;
}

// Hyperbolic tangent function
template <int Order>
Dual<Order> tanh(const Dual<Order>& x) {
  Dual<Order> result(std::tanh(x.val), x.der.size(), -1);
    const double sech2_val = 1.0 - result.val * result.val;
    for (size_t i = 0; i < x.der.size(); ++i) {
        result.der[i] = (1 - sech2_val) * x.der[i];
    }
    return result;
}

#endif
