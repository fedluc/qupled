#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include "num_util.hpp"
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <memory>
#include <vector>

// -----------------------------------------------------------------
// C++ wrappers to GSL objects
// -----------------------------------------------------------------

namespace GslWrappers {

  // Wrapper to gsl_function
  template <typename T>
  class GslFunctionWrap : public gsl_function {

  private:

    const T &func;
    static double invoke(double x, void *params) {
      return static_cast<GslFunctionWrap *>(params)->func(x);
    }

  public:

    explicit GslFunctionWrap(const T &func_)
        : func(func_) {
      function = &GslFunctionWrap::invoke;
      params = this;
    }
  };

  // Wrapper to handle errors in GSL functions
  template <typename Func, typename... Args>
  void callGSLFunction(Func &&gslFunction, Args &&...args);

  // Wrapper to handle allocation errors in GSL functions
  template <typename Ptr, typename Func, typename... Args>
  void callGSLAlloc(Ptr &ptr, Func &&gslFunction, Args &&...args);

} // namespace GslWrappers

// -----------------------------------------------------------------
// Classes to interpolate data
// -----------------------------------------------------------------

// Interpolator for 1D data
class Interpolator1D {

public:

  // Constructor
  Interpolator1D(const std::vector<double> &x, const std::vector<double> &y);
  Interpolator1D(const double &x, const double &y, const size_t n_);
  explicit Interpolator1D();
  // Destructor
  ~Interpolator1D();
  // Check if the data used to construct the interpolator is valid
  bool isValid() const;
  // Reset
  void reset(const double &x, const double &y, const size_t n_);
  // Evaluate
  double eval(const double &x) const;

private:

  // Type
  const gsl_interp_type *TYPE = gsl_interp_cspline;
  // Spline
  gsl_spline *spline;
  // Accelerator
  gsl_interp_accel *acc;
  // Cutoff (extrapolation for x > cutoff)
  double cutoff;
  // Size
  size_t n;
  //
  // Setup interpolator
  void setup(const double &x, const double &y, const size_t n_);
};

// Interpolator for 2D data
class Interpolator2D {

public:

  // Constructor
  Interpolator2D(const double &x,
                 const double &y,
                 const double &z,
                 const int nx_,
                 const int ny_);
  explicit Interpolator2D(const Interpolator2D &it);
  explicit Interpolator2D();
  // Destructor
  ~Interpolator2D();
  // Check if the data used to construct the interpolator is valid
  bool isValid() const;
  // Reset
  void reset(const double &x,
             const double &y,
             const double &z,
             const int szx_,
             const int szy_);
  // Evaluate
  double eval(const double &x, const double &y) const;

private:

  // Type
  const gsl_interp2d_type *TYPE = gsl_interp2d_bicubic;
  // Spline
  gsl_spline2d *spline;
  // Accelerator
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
  // Size
  size_t nx;
  size_t ny;
  // Setup interpolator
  void setup(const double &x,
             const double &y,
             const double &z,
             const int nx_,
             const int ny_);
};

// -----------------------------------------------------------------
// Classes to find roots of equations
// -----------------------------------------------------------------

class RootSolverBase {

public:

  double getSolution() const { return sol; };

protected:

  // Accuracy
  const double relErr;
  // Iterations
  const int maxIter;
  int iter;
  // Solver status
  int status;
  // Solution
  double sol;
  // Protected constructor
  RootSolverBase(const double &relErr_, const int maxIter_)
      : relErr(relErr_),
        maxIter(maxIter_),
        iter(0),
        status(GSL_CONTINUE) {}
  explicit RootSolverBase()
      : RootSolverBase(1.0e-10, 1000) {}
};

class BrentRootSolver : public RootSolverBase {

public:

  explicit BrentRootSolver();
  ~BrentRootSolver();
  void solve(const std::function<double(double)> &func,
             const std::vector<double> &guess);

private:

  // Function to solve
  gsl_function *F;
  // Type of solver
  const gsl_root_fsolver_type *rst;
  // Solver
  gsl_root_fsolver *rs;
};

class SecantSolver : public RootSolverBase {

public:

  SecantSolver(const double relErr_, const int maxIter_)
      : RootSolverBase(relErr_, maxIter_) {}
  explicit SecantSolver() {}
  void solve(const std::function<double(double)> &func,
             const std::vector<double> &guess);
};

// -----------------------------------------------------------------
// Class to compute 1D integrals
// -----------------------------------------------------------------

class Integrator1D {

public:

  // Typdef
  enum Type { DEFAULT, FOURIER, SINGULAR };

  class Param {

  public:

    const double xMin = numUtil::NaN;
    const double xMax = numUtil::NaN;
    const double fourierR = numUtil::NaN;
    Param(const double &xMin_, const double &xMax_)
        : xMin(xMin_),
          xMax(xMax_) {}
    explicit Param(const double &fourierR_)
        : fourierR(fourierR_) {}
  };

  // Constructors
  Integrator1D(const Type &type, const double &relErr);
  explicit Integrator1D(const double &relErr)
      : Integrator1D(Type::DEFAULT, relErr) {}
  Integrator1D(const Integrator1D &other)
      : Integrator1D(other.getType(), other.getAccuracy()) {}
  // Compute integral
  void compute(const std::function<double(double)> &func,
               const Param &param) const;
  // Getters
  double getSolution() const;
  double getAccuracy() const { return gslIntegrator->getAccuracy(); }
  Type getType() const { return gslIntegrator->getType(); }

private:

  // Base class for all integrators derived from GSL
  class Base {

  public:

    // Constructors
    Base(const Type &type_, const size_t &limit_, const double &relErr_)
        : type(type_),
          limit(limit_),
          relErr(relErr_) {}
    // Destructor
    virtual ~Base() = default;
    // Getters
    double getSolution() const { return sol; }
    double getAccuracy() const { return relErr; }
    Type getType() const { return type; }
    // Compute integral
    virtual void compute(const std::function<double(double)> &func,
                         const Param &param) = 0;

  protected:

    // Function to integrate
    gsl_function *F;
    // Integrator type
    const Type type;
    // Integration workspace limit
    const size_t limit;
    // Accuracy
    const double relErr;
    // Residual error
    double err;
    // Solution
    double sol;
  };

  // CQUAD integrator from GSL
  class CQUAD : public Base {
  public:

    // Constructors
    explicit CQUAD(const double &relErr_);
    CQUAD(const CQUAD &other)
        : Integrator1D::CQUAD(other.relErr) {}
    // Destructor
    ~CQUAD();
    // Compute integral
    void compute(const std::function<double(double)> &func,
                 const Param &param) override;

  private:

    // Integration workspace
    gsl_integration_cquad_workspace *wsp;
    // Number of evaluations
    size_t nEvals;
  };

  // QAWO integrator from GSL
  class QAWO : public Base {

  public:

    // Constructors
    explicit QAWO(const double &relErr_);
    QAWO(const QAWO &other)
        : Integrator1D::QAWO(other.relErr) {}
    // Destructor
    ~QAWO();
    // Compute integral
    void compute(const std::function<double(double)> &func,
                 const Param &param) override;

  private:

    // Integration workspace
    gsl_integration_workspace *wsp;
    gsl_integration_workspace *wspc;
    gsl_integration_qawo_table *qtab;
  };

  // QAGS integrator from GSL
  class QAGS : public Base {

  public:

    // Constructors
    explicit QAGS(const double &relErr_);
    QAGS(const QAGS &other)
        : Integrator1D::QAGS(other.relErr) {}
    // Destructor
    ~QAGS();
    // Compute integral
    void compute(const std::function<double(double)> &func,
                 const Param &param) override;

  private:

    // Integration workspace
    gsl_integration_workspace *wsp;
  };

  // Pointers to GSL integrals
  std::unique_ptr<Base> gslIntegrator;
};

// -----------------------------------------------------------------
// Class to compute 2D integrals
// -----------------------------------------------------------------

class Integrator2D {

public:

  // Typedef
  using Type = Integrator1D::Type;
  using Param1D = Integrator1D::Param;
  // Class to handle integration parameters
  class Param : public Param1D {

  public:

    using Func = std::function<double(double)>;
    const Func yMin = [&](const double &x) {
      (void)(x);
      return yMinNum;
    };
    const Func yMax = [&](const double &x) {
      (void)(x);
      return yMaxNum;
    };
    const std::vector<double> xGrid;
    Param(const double &xMin_,
          const double &xMax_,
          const Func &yMin_,
          const Func &yMax_)
        : Integrator1D::Param(xMin_, xMax_),
          yMin(yMin_),
          yMax(yMax_) {}
    Param(const double &xMin_,
          const double &xMax_,
          const double &yMin_,
          const double &yMax_)
        : Integrator1D::Param(xMin_, xMax_),
          yMinNum(yMin_),
          yMaxNum(yMax_) {}
    Param(const double &fourierR_)
        : Integrator1D::Param(fourierR_) {}

  private:

    const double yMinNum = numUtil::NaN;
    const double yMaxNum = numUtil::NaN;
  };

  // Constructors
  Integrator2D(const Type &type1, const Type &type2, const double &relErr)
      : itg1(type1, relErr),
        itg2(type2, relErr) {}
  Integrator2D(const Type &type, const double &relErr)
      : Integrator2D(type, type, relErr) {}
  Integrator2D(const double &relErr)
      : Integrator2D(Type::DEFAULT, relErr) {}
  // Compute integral
  void compute(const std::function<double(double)> &func1,
               const std::function<double(double)> &func2,
               const Param &param,
               const std::vector<double> &xGrid);
  // Getters
  double getX() const { return x; };
  double getSolution() const { return sol; };

private:

  // Level 1 integrator (outermost integral)
  Integrator1D itg1;
  // Level 2 integrator
  Integrator1D itg2;
  // Temporary variable for level 2 integration
  double x;
  // Solution
  double sol;
};

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

class Dual2 {
public:

  double val;
  double dx;
  double dy;
  double dxx;
  double dyy;
  double dxy;

  // Constructor
  Dual2(double value,
        double dx = 0.0,
        double dy = 0.0,
        double dxx = 0.0,
        double dyy = 0.0,
        double dxy = 0.0)
      : val(value),
        dx(dx),
        dy(dy),
        dxx(dxx),
        dyy(dyy),
        dxy(dxy) {}

  // Overload addition operator
  Dual2 operator+(const Dual2 &other) const {
    return Dual2(val + other.val,
                 dx + other.dx,
                 dy + other.dy,
                 dxx + other.dxx,
                 dyy + other.dyy,
                 dxy + other.dxy);
  }

  // Overload subtraction operator
  Dual2 operator-(const Dual2 &other) const {
    return Dual2(val - other.val,
                 dx - other.dx,
                 dy - other.dy,
                 dxx - other.dxx,
                 dyy - other.dyy,
                 dxy - other.dxy);
  }

  // Overload multiplication operator
  Dual2 operator*(const Dual2 &other) const {
    return Dual2(val * other.val,
                 val * other.dx + dx * other.val,
                 val * other.dy + dy * other.val,
                 val * other.dxx + 2 * dx * other.dx + dxx * other.val,
                 val * other.dyy + 2 * dy * other.dy + dyy * other.val,
                 val * other.dxy + dx * other.dy + dy * other.dx +
                     dxy * other.val);
  }

  // Overload division operator
  Dual2 operator/(const Dual2 &other) const {
    const double inv_val = 1.0 / other.val;
    const double inv_val2 = inv_val * inv_val;
    return Dual2(
        val * inv_val,
        (dx - val * other.dx * inv_val) * inv_val,
        (dy - val * other.dy * inv_val) * inv_val,
        (dxx - 2 * dx * other.dx * inv_val - val * other.dxx * inv_val +
         2 * val * other.dx * other.dx * inv_val2) *
            inv_val,
        (dyy - 2 * dy * other.dy * inv_val - val * other.dyy * inv_val +
         2 * val * other.dy * other.dy * inv_val2) *
            inv_val,
        (dxy - dx * other.dy * inv_val - dy * other.dx * inv_val -
         val * other.dxy * inv_val + 2 * val * other.dx * other.dy * inv_val2) *
            inv_val);
  }

  // Overload addition operator: Dual2 + double
  Dual2 operator+(double scalar) const {
    return Dual2(val + scalar, dx, dy, dxx, dyy, dxy);
  }

  // Overload addition operator: double + Dual2
  friend Dual2 operator+(double scalar, const Dual2 &dual) {
    return dual + scalar; // Reuse the Dual2 + double operator
  }

  // Overload subtraction operator: Dual2 - double
  Dual2 operator-(double scalar) const {
    return Dual2(val - scalar, dx, dy, dxx, dyy, dxy);
  }

  // Overload subtraction operator: double - Dual2
  friend Dual2 operator-(double scalar, const Dual2 &dual) {
    return Dual2(
        scalar - dual.val, -dual.dx, -dual.dy, -dual.dxx, -dual.dyy, -dual.dxy);
  }

  // Overload multiplication operator: Dual2 * double
  Dual2 operator*(double scalar) const {
    return Dual2(val * scalar,
                 dx * scalar,
                 dy * scalar,
                 dxx * scalar,
                 dyy * scalar,
                 dxy * scalar);
  }

  // Overload multiplication operator: double * Dual2
  friend Dual2 operator*(double scalar, const Dual2 &dual) {
    return dual * scalar; // Reuse the Dual2 * double operator
  }

  // Overload division operator: Dual2 / double
  Dual2 operator/(double scalar) const {
    const double inv_scalar = 1.0 / scalar;
    return Dual2(val * inv_scalar,
                 dx * inv_scalar,
                 dy * inv_scalar,
                 dxx * inv_scalar,
                 dyy * inv_scalar,
                 dxy * inv_scalar);
  }

  // Overload division operator: double / Dual2
  friend Dual2 operator/(double scalar, const Dual2 &dual) {
    const double inv_val = 1.0 / dual.val;
    return Dual2(scalar * inv_val,
                 -scalar * dual.dx * inv_val * inv_val,
                 -scalar * dual.dy * inv_val * inv_val,
                 scalar * (2 * dual.dx * dual.dx - dual.dxx * dual.val) *
                     inv_val * inv_val * inv_val,
                 scalar * (2 * dual.dy * dual.dy - dual.dyy * dual.val) *
                     inv_val * inv_val * inv_val,
                 scalar * (2 * dual.dx * dual.dy - dual.dxy * dual.val) *
                     inv_val * inv_val * inv_val);
  }

  // Overload sine function
  friend Dual2 sin(const Dual2 &x) {
    return Dual2(std::sin(x.val),
                 std::cos(x.val) * x.dx,
                 std::cos(x.val) * x.dy,
                 -std::sin(x.val) * x.dx * x.dx + std::cos(x.val) * x.dxx -
                     std::sin(x.val) * x.dy * x.dy + std::cos(x.val) * x.dyy -
                     std::sin(x.val) * x.dx * x.dy + std::cos(x.val) * x.dxy);
  }

  // Overload exponential function
  friend Dual2 exp(const Dual2 &x) {
    const double exp_val = std::exp(x.val);
    return Dual2(exp_val,
                 exp_val * x.dx,
                 exp_val * x.dy,
                 exp_val * (x.dx * x.dx + x.dxx),
                 exp_val * (x.dy * x.dy + x.dyy),
                 exp_val * (x.dx * x.dy + x.dxy));
  }

  // Overload square root function
  friend Dual2 sqrt(const Dual2 &x) {
    const double sqrt_val = std::sqrt(x.val);
    const double inv_sqrt = 0.5 / sqrt_val;
    return Dual2(sqrt_val,
                 x.dx * inv_sqrt,
                 x.dy * inv_sqrt,
                 (x.dxx - x.dx * x.dx / (2 * x.val)) * inv_sqrt,
                 (x.dyy - x.dy * x.dy / (2 * x.val)) * inv_sqrt,
                 (x.dxy - x.dx * x.dy / (2 * x.val)) * inv_sqrt);
  }

  // Overload hyperbolic tangent function
  friend Dual2 tanh(const Dual2 &x) {
    const double tanh_val = std::tanh(x.val);
    const double sech2_val = 1.0 - tanh_val * tanh_val;
    return Dual2(tanh_val,
                 x.dx * sech2_val,
                 x.dy * sech2_val,
                 sech2_val * (x.dxx - 2 * x.dx * tanh_val * x.dx),
                 sech2_val * (x.dyy - 2 * x.dy * tanh_val * x.dy),
                 sech2_val * (x.dxy - tanh_val * (x.dx * x.dy + x.dy * x.dx)));
  }
};

#endif
