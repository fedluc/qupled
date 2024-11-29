#include <cmath>
#include "auto_diff.hpp"

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

// addition operators
AutoDiff2 operator+(const AutoDiff2 &dual1, const AutoDiff2 &dual2) {
    return AutoDiff2(dual1.val + dual2.val,
                 dual1.dx + dual2.dx,
                 dual1.dy + dual2.dy,
                 dual1.dxx + dual2.dxx,
                 dual1.dyy + dual2.dyy,
                 dual1.dxy + dual2.dxy);
}

AutoDiff2 operator+(const AutoDiff2 &dual, const double &scalar) {
    return AutoDiff2(dual.val + scalar, dual.dx, dual.dy, dual.dxx, dual.dyy, dual.dxy);
}

AutoDiff2 operator+(const double &scalar, const AutoDiff2 &dual) {
    return dual + scalar;
}

// subtraction operators
AutoDiff2 operator-(const AutoDiff2 &dual1, const AutoDiff2 &dual2) {
    return AutoDiff2(dual1.val - dual2.val,
                 dual1.dx - dual2.dx,
                 dual1.dy - dual2.dy,
                 dual1.dxx - dual2.dxx,
                 dual1.dyy - dual2.dyy,
                 dual1.dxy - dual2.dxy);
}

AutoDiff2 operator-(const AutoDiff2 &dual, double scalar) {
    return AutoDiff2(dual.val - scalar, dual.dx, dual.dy, dual.dxx, dual.dyy, dual.dxy);
}

AutoDiff2 operator-(const double &scalar, const AutoDiff2 &dual) {
    return AutoDiff2(scalar - dual.val, -dual.dx, -dual.dy, -dual.dxx, -dual.dyy, -dual.dxy);
}

// multiplication operators
AutoDiff2 operator*(const AutoDiff2 &dual1, const AutoDiff2 &dual2) {
    return AutoDiff2(dual1.val * dual2.val,
                 dual1.val * dual2.dx + dual1.dx * dual2.val,
                 dual1.val * dual2.dy + dual1.dy * dual2.val,
                 dual1.val * dual2.dxx + 2 * dual1.dx * dual2.dx + dual1.dxx * dual2.val,
                 dual1.val * dual2.dyy + 2 * dual1.dy * dual2.dy + dual1.dyy * dual2.val,
                 dual1.val * dual2.dxy + dual1.dx * dual2.dy + dual1.dy * dual2.dx + dual1.dxy * dual2.val);
}

AutoDiff2 operator*(const AutoDiff2 &dual, const double &scalar) {
    return AutoDiff2(dual.val * scalar, dual.dx * scalar, dual.dy * scalar, dual.dxx * scalar, dual.dyy * scalar, dual.dxy * scalar);
}

AutoDiff2 operator*(const double &scalar, const AutoDiff2 &dual) {
    return dual * scalar;
}

// division operators
AutoDiff2 operator/(const AutoDiff2 &dual1, const AutoDiff2 &dual2) {
    const double inv_val = 1.0 / dual2.val;
    const double inv_val2 = inv_val * inv_val;
    return AutoDiff2(
        dual1.val * inv_val,
        (dual1.dx - dual1.val * dual2.dx * inv_val) * inv_val,
        (dual1.dy - dual1.val * dual2.dy * inv_val) * inv_val,
        (dual1.dxx - 2 * dual1.dx * dual2.dx * inv_val - dual1.val * dual2.dxx * inv_val +
         2 * dual1.val * dual2.dx * dual2.dx * inv_val2) *
            inv_val,
        (dual1.dyy - 2 * dual1.dy * dual2.dy * inv_val - dual1.val * dual2.dyy * inv_val +
         2 * dual1.val * dual2.dy * dual2.dy * inv_val2) *
            inv_val,
        (dual1.dxy - dual1.dx * dual2.dy * inv_val - dual1.dy * dual2.dx * inv_val -
         dual1.val * dual2.dxy * inv_val + 2 * dual1.val * dual2.dx * dual2.dy * inv_val2) *
            inv_val);
}

AutoDiff2 operator/(const AutoDiff2 &dual, const double &scalar) {
    const double inv_scalar = 1.0 / scalar;
    return AutoDiff2(dual.val * inv_scalar,
                 dual.dx * inv_scalar,
                 dual.dy * inv_scalar,
                 dual.dxx * inv_scalar,
                 dual.dyy * inv_scalar,
                 dual.dxy * inv_scalar);
}

AutoDiff2 operator/(const double &scalar, const AutoDiff2 &dual) {
    const double inv_val = 1.0 / dual.val;
    return AutoDiff2(scalar * inv_val,
                 -scalar * dual.dx * inv_val * inv_val,
                 -scalar * dual.dy * inv_val * inv_val,
                 scalar * (2 * dual.dx * dual.dx - dual.dxx * dual.val) * inv_val * inv_val * inv_val,
                 scalar * (2 * dual.dy * dual.dy - dual.dyy * dual.val) * inv_val * inv_val * inv_val,
                 scalar * (2 * dual.dx * dual.dy - dual.dxy * dual.val) * inv_val * inv_val * inv_val);
}

// exponential function
AutoDiff2 exp(const AutoDiff2 &x) {
    const double exp_val = std::exp(x.val);
    return AutoDiff2(exp_val,
                 exp_val * x.dx,
                 exp_val * x.dy,
                 exp_val * (x.dx * x.dx + x.dxx),
                 exp_val * (x.dy * x.dy + x.dyy),
                 exp_val * (x.dx * x.dy + x.dxy));
}

// square root function
AutoDiff2 sqrt(const AutoDiff2 &x) {
    const double sqrt_val = std::sqrt(x.val);
    const double inv_sqrt = 0.5 / sqrt_val;
    return AutoDiff2(sqrt_val,
                 x.dx * inv_sqrt,
                 x.dy * inv_sqrt,
                 (x.dxx - x.dx * x.dx / (2 * x.val)) * inv_sqrt,
                 (x.dyy - x.dy * x.dy / (2 * x.val)) * inv_sqrt,
                 (x.dxy - x.dx * x.dy / (2 * x.val)) * inv_sqrt);
}

// hyperbolic tangent function
AutoDiff2 tanh(const AutoDiff2 &x) {
    const double tanh_val = std::tanh(x.val);
    const double sech2_val = 1.0 - tanh_val * tanh_val;
    return AutoDiff2(tanh_val,
                 x.dx * sech2_val,
                 x.dy * sech2_val,
                 sech2_val * (x.dxx - 2 * x.dx * tanh_val * x.dx),
                 sech2_val * (x.dyy - 2 * x.dy * tanh_val * x.dy),
                 sech2_val * (x.dxy - tanh_val * (x.dx * x.dy + x.dy * x.dx)));
}
