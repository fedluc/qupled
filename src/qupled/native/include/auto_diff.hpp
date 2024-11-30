#ifndef AUTO_DIFF_HPP
#define AUTO_DIFF_HPP

#include <cmath>
#include <vector>

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

// Automatic differentiation for two dimensional functions
class AutoDiff1 {
public:

  double val;
  double dx;
  double dy;

  // Constructors
  AutoDiff1(double val_, double dx_ = 0.0, double dy_ = 0.0)
      : val(val_),
        dx(dx_),
        dy(dy_) {}
};

// addition operators
AutoDiff1 operator+(const AutoDiff1 &dual1, const AutoDiff1 &dual2);
AutoDiff1 operator+(const AutoDiff1 &dual, const double &scalar);
AutoDiff1 operator+(const double &scalar, const AutoDiff1 &dual);

// subtraction operators
AutoDiff1 operator-(const AutoDiff1 &dual1, const AutoDiff1 &dual2);
AutoDiff1 operator-(const AutoDiff1 &dual, double scalar);
AutoDiff1 operator-(const double &scalar, const AutoDiff1 &dual);

// multiplication operators
AutoDiff1 operator*(const AutoDiff1 &dual1, const AutoDiff1 &dual2);
AutoDiff1 operator*(const AutoDiff1 &dual, const double &scalar);
AutoDiff1 operator*(const double &scalar, const AutoDiff1 &dual);

// division operators
AutoDiff1 operator/(const AutoDiff1 &dual1, const AutoDiff1 &dual2);
AutoDiff1 operator/(const AutoDiff1 &dual, const double &scalar);
AutoDiff1 operator/(const double &scalar, const AutoDiff1 &dual);

// exponential function
AutoDiff1 exp(const AutoDiff1 &x);

// square root function
AutoDiff1 sqrt(const AutoDiff1 &x);

// hyperbolic tangent function
AutoDiff1 tanh(const AutoDiff1 &x);

// Automatic differentiation for two dimensional functions
class AutoDiff2 {
public:

  AutoDiff1 val;
  AutoDiff1 dx;
  AutoDiff1 dy;

  AutoDiff2(const AutoDiff1 &val_, const AutoDiff1 &dx_, const AutoDiff1 &dy_)
      : val(val_),
        dx(dx_),
        dy(dy_) {}

  AutoDiff2(const double &val_,
            const double &dx_ = 0.0,
            const double &dy_ = 0.0)
      : AutoDiff2(AutoDiff1(val_, dx_, dy_), AutoDiff1(dx_), AutoDiff1(dy_)) {}
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

#endif
