#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include "util/num_util.hpp"
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <memory>
#include <vector>

/**
 * @brief Thin C++ wrappers around GSL objects.
 *
 * Provides type-safe adapters for @c gsl_function (function objects) and
 * error-checked wrappers for GSL allocation and call patterns.
 */
namespace GslWrappers {

  /**
   * @brief Adapts a callable object to a @c gsl_function.
   *
   * @tparam T Type of the callable (lambda, functor, or function pointer).
   */
  template <typename T>
  class GslFunctionWrap : public gsl_function {

  private:

    /** @brief The wrapped callable. */
    const T &func;

    /**
     * @brief Static trampoline invoked by GSL.
     * @param x      Evaluation point.
     * @param params Pointer to the @p GslFunctionWrap instance.
     * @return @f$f(x)@f$.
     */
    static double invoke(double x, void *params) {
      return static_cast<GslFunctionWrap *>(params)->func(x);
    }

  public:

    /**
     * @brief Construct the wrapper.
     * @param func_ Callable to wrap.
     */
    explicit GslFunctionWrap(const T &func_)
        : func(func_) {
      function = &GslFunctionWrap::invoke;
      params = this;
    }
  };

  /**
   * @brief Call a GSL function and propagate any GSL error as a C++ exception.
   *
   * @tparam Func    GSL function type.
   * @tparam Args    Argument types.
   * @param gslFunction GSL function to call.
   * @param args        Arguments forwarded to @p gslFunction.
   */
  template <typename Func, typename... Args>
  void callGSLFunction(Func &&gslFunction, Args &&...args);

  /**
   * @brief Allocate a GSL resource and propagate any allocation error.
   *
   * @tparam Ptr     Pointer type to receive the allocation result.
   * @tparam Func    GSL allocation function type.
   * @tparam Args    Argument types.
   * @param ptr         Pointer to populate with the allocated resource.
   * @param gslFunction GSL allocation function.
   * @param args        Arguments forwarded to @p gslFunction.
   */
  template <typename Ptr, typename Func, typename... Args>
  void callGSLAlloc(Ptr &ptr, Func &&gslFunction, Args &&...args);

} // namespace GslWrappers

/**
 * @brief Wrappers for commonly needed GSL special functions.
 */
namespace SpecialFunctions {

  /**
   * @brief Fermi–Dirac integral of order 1/2.
   * @param x Argument.
   * @return @f$F_{1/2}(x) = \frac{1}{\Gamma(3/2)} \int_0^\infty
   * \frac{t^{1/2}}{e^{t-x} + 1} dt@f$.
   */
  double fermiDirac12(const double &x);

  /**
   * @brief Fermi–Dirac integral of order -1/2.
   * @param x Argument.
   * @return @f$F_{-1/2}(x) = \frac{1}{\Gamma(1/2)} \int_0^\infty
   * \frac{t^{-1/2}}{e^{t-x} + 1} dt@f$.
   */
  double fermiDiracm12(const double &x);

  /**
   * @brief Hyperbolic cotangent.
   * @param x Argument.
   * @return @f$\coth(x)@f$.
   */
  double coth(const double &x);

  /**
   * @brief Complete elliptic integral of the first kind.
   * @param x Modulus @f$k@f$.
   * @return @f$K(k)@f$.
   */
  double ellipticK(const double &x);

  /**
   * @brief Complete elliptic integral of the second kind.
   * @param x Modulus @f$k@f$.
   * @return @f$E(k)@f$.
   */
  double ellipticE(const double &x);

  /**
   * @brief Bessel function of the first kind of order zero.
   * @param x Argument.
   * @return @f$J_0(x)@f$.
   */
  double besselJ0(const double &x);

} // namespace SpecialFunctions

/**
 * @brief 1D cubic-spline interpolator backed by GSL.
 *
 * Wraps @c gsl_spline with caching of the last evaluation point via a
 * @c gsl_interp_accel. Values outside the data range are extrapolated
 * by holding the boundary value (constant extrapolation beyond @p cutoff).
 */
class Interpolator1D {

public:

  /**
   * @brief Construct from separate x and y data vectors.
   * @param x Strictly increasing abscissa values.
   * @param y Ordinate values (same length as @p x).
   */
  Interpolator1D(const std::vector<double> &x, const std::vector<double> &y);

  /**
   * @brief Construct from raw pointers.
   * @param x  Pointer to the first abscissa value.
   * @param y  Pointer to the first ordinate value.
   * @param n_ Number of data points.
   */
  Interpolator1D(const double &x, const double &y, const size_t n_);

  /** @brief Construct an empty (invalid) interpolator. */
  explicit Interpolator1D();

  /** @brief Destructor — releases GSL resources. */
  ~Interpolator1D();

  /**
   * @brief Return true if the interpolator holds valid data.
   * @return False for default-constructed or reset-but-not-populated objects.
   */
  bool isValid() const;

  /**
   * @brief Reinitialize with new raw-pointer data.
   * @param x  Pointer to the first abscissa value.
   * @param y  Pointer to the first ordinate value.
   * @param n_ Number of data points.
   */
  void reset(const double &x, const double &y, const size_t n_);

  /**
   * @brief Evaluate the spline at @p x.
   * @param x Evaluation point.
   * @return Interpolated value (constant extrapolation beyond the data range).
   */
  double eval(const double &x) const;

private:

  /** @brief GSL spline type (cubic). */
  const gsl_interp_type *TYPE = gsl_interp_cspline;
  /** @brief GSL spline object. */
  gsl_spline *spline;
  /** @brief GSL interpolation accelerator. */
  gsl_interp_accel *acc;
  /** @brief Maximum abscissa value of the data range. */
  double cutoff;
  /** @brief Number of data points. */
  size_t n;
  /**
   * @brief Shared initialization helper.
   * @param x  Pointer to the first abscissa value.
   * @param y  Pointer to the first ordinate value.
   * @param n_ Number of data points.
   */
  void setup(const double &x, const double &y, const size_t n_);
};

/**
 * @brief 2D bicubic-spline interpolator backed by GSL.
 *
 * Wraps @c gsl_spline2d. The data must be on a regular rectangular grid
 * with @p nx abscissa values and @p ny ordinate values.
 */
class Interpolator2D {

public:

  /**
   * @brief Construct from raw-pointer grid data.
   * @param x   Pointer to the first x abscissa value.
   * @param y   Pointer to the first y abscissa value.
   * @param z   Pointer to the first grid value (row-major: z[i*ny + j]).
   * @param nx_ Number of x grid points.
   * @param ny_ Number of y grid points.
   */
  Interpolator2D(const double &x,
                 const double &y,
                 const double &z,
                 const int nx_,
                 const int ny_);

  /**
   * @brief Copy constructor (deep copy of the underlying GSL objects).
   * @param it Source interpolator.
   */
  explicit Interpolator2D(const Interpolator2D &it);

  /** @brief Construct an empty (invalid) interpolator. */
  explicit Interpolator2D();

  /** @brief Destructor — releases GSL resources. */
  ~Interpolator2D();

  /**
   * @brief Return true if the interpolator holds valid data.
   * @return False for default-constructed objects.
   */
  bool isValid() const;

  /**
   * @brief Reinitialize with new grid data.
   * @param x    Pointer to the first x abscissa.
   * @param y    Pointer to the first y abscissa.
   * @param z    Pointer to the first grid value.
   * @param szx_ Number of x grid points.
   * @param szy_ Number of y grid points.
   */
  void reset(const double &x,
             const double &y,
             const double &z,
             const int szx_,
             const int szy_);

  /**
   * @brief Evaluate the bicubic spline at (@p x, @p y).
   * @param x x-coordinate.
   * @param y y-coordinate.
   * @return Interpolated value.
   */
  double eval(const double &x, const double &y) const;

private:

  /** @brief GSL 2D spline type (bicubic). */
  const gsl_interp2d_type *TYPE = gsl_interp2d_bicubic;
  /** @brief GSL 2D spline object. */
  gsl_spline2d *spline;
  /** @brief GSL accelerator for the x axis. */
  gsl_interp_accel *xacc;
  /** @brief GSL accelerator for the y axis. */
  gsl_interp_accel *yacc;
  /** @brief Number of x grid points. */
  size_t nx;
  /** @brief Number of y grid points. */
  size_t ny;
  /**
   * @brief Shared initialization helper.
   * @param x    Pointer to the first x abscissa.
   * @param y    Pointer to the first y abscissa.
   * @param z    Pointer to the first grid value.
   * @param nx_  Number of x grid points.
   * @param ny_  Number of y grid points.
   */
  void setup(const double &x,
             const double &y,
             const double &z,
             const int nx_,
             const int ny_);
};

/**
 * @brief Base class for scalar root solvers.
 *
 * Tracks iteration count, solver status, and the solution in a uniform way
 * for all concrete root-solver implementations.
 */
class RootSolverBase {

public:

  /**
   * @brief Return the last computed solution.
   * @return Root value found by the solver.
   */
  double getSolution() const { return sol; };

protected:

  /** @brief Relative error tolerance for convergence. */
  const double relErr;
  /** @brief Maximum number of iterations. */
  const int maxIter;
  /** @brief Current iteration count. */
  int iter;
  /** @brief GSL solver status code (@c GSL_CONTINUE until converged). */
  int status;
  /** @brief Root found by the solver. */
  double sol;

  /**
   * @brief Construct with explicit tolerances.
   * @param relErr_   Relative error tolerance.
   * @param maxIter_  Maximum iteration count.
   */
  RootSolverBase(const double &relErr_, const int maxIter_)
      : relErr(relErr_),
        maxIter(maxIter_),
        iter(0),
        status(GSL_CONTINUE) {}

  /** @brief Construct with default tolerances (1e-10, 1000). */
  explicit RootSolverBase()
      : RootSolverBase(1.0e-10, 1000) {}
};

/**
 * @brief Root solver using the GSL Brent bracketing algorithm.
 *
 * Requires an initial bracket @f$[a, b]@f$ with @f$f(a) \cdot f(b) < 0@f$.
 */
class BrentRootSolver : public RootSolverBase {

public:

  /** @brief Construct the Brent solver with default tolerances. */
  explicit BrentRootSolver();

  /** @brief Destructor — releases GSL resources. */
  ~BrentRootSolver();

  /**
   * @brief Find the root of @p func within the initial bracket @p guess.
   * @param func  Function whose root is sought.
   * @param guess Two-element vector @f$[a, b]@f$ bracketing the root.
   */
  void solve(const std::function<double(double)> &func,
             const std::vector<double> &guess);

private:

  /** @brief GSL function wrapper. */
  gsl_function *F;
  /** @brief GSL root-solver algorithm type. */
  const gsl_root_fsolver_type *rst;
  /** @brief GSL root-solver workspace. */
  gsl_root_fsolver *rs;
};

/**
 * @brief Root solver using the secant (chord) method.
 *
 * Does not require a bracket but needs two initial guesses that are
 * sufficiently close to the root.
 */
class SecantSolver : public RootSolverBase {

public:

  /**
   * @brief Construct with explicit tolerances.
   * @param relErr_  Relative error tolerance.
   * @param maxIter_ Maximum iteration count.
   */
  SecantSolver(const double relErr_, const int maxIter_)
      : RootSolverBase(relErr_, maxIter_) {}

  /** @brief Construct with default tolerances. */
  explicit SecantSolver() {}

  /**
   * @brief Find the root of @p func starting from the two guesses in @p guess.
   * @param func  Function whose root is sought.
   * @param guess Two-element vector of initial approximations.
   */
  void solve(const std::function<double(double)> &func,
             const std::vector<double> &guess);
};

/**
 * @brief 1D numerical integrator wrapping GSL quadrature routines.
 *
 * Supports three quadrature modes (CQUAD for general integrands, QAWO for
 * oscillatory Fourier-type integrals, QAGS for integrands with endpoint
 * singularities) selected via the @p Type enum at construction time.
 */
class Integrator1D {

public:

  /**
   * @brief Quadrature mode selection.
   *
   * - @p DEFAULT: GSL CQUAD (general-purpose doubly adaptive).
   * - @p FOURIER: GSL QAWO (weighted oscillatory integrals).
   * - @p SINGULAR: GSL QAGS (endpoint-singularity handling).
   */
  enum Type { DEFAULT, FOURIER, SINGULAR };

  /**
   * @brief Integration parameters.
   *
   * Encapsulates the integration limits for finite intervals or the Fourier
   * radius for QAWO-type integrals.
   */
  class Param {

  public:

    /** @brief Lower integration limit (NaN for FOURIER type). */
    const double xMin = numUtil::NaN;
    /** @brief Upper integration limit (NaN for FOURIER type). */
    const double xMax = numUtil::NaN;
    /** @brief Fourier integration radius (NaN for finite-interval types). */
    const double fourierR = numUtil::NaN;

    /**
     * @brief Construct for a finite-interval integral.
     * @param xMin_ Lower limit.
     * @param xMax_ Upper limit.
     */
    Param(const double &xMin_, const double &xMax_)
        : xMin(xMin_),
          xMax(xMax_) {}

    /**
     * @brief Construct for a Fourier-type integral.
     * @param fourierR_ Fourier radius.
     */
    explicit Param(const double &fourierR_)
        : fourierR(fourierR_) {}
  };

  /**
   * @brief Construct with an explicit quadrature type and relative error.
   * @param type    Quadrature mode.
   * @param relErr  Relative error tolerance.
   */
  Integrator1D(const Type &type, const double &relErr);

  /**
   * @brief Construct with the default (CQUAD) quadrature type.
   * @param relErr Relative error tolerance.
   */
  explicit Integrator1D(const double &relErr)
      : Integrator1D(Type::DEFAULT, relErr) {}

  /**
   * @brief Copy constructor.
   * @param other Source integrator; type and accuracy are copied.
   */
  Integrator1D(const Integrator1D &other)
      : Integrator1D(other.getType(), other.getAccuracy()) {}

  /**
   * @brief Evaluate the integral of @p func over the domain specified by @p
   * param.
   * @param func  Integrand function @f$f : \mathbb{R} \to \mathbb{R}@f$.
   * @param param Integration parameters (limits or Fourier radius).
   */
  void compute(const std::function<double(double)> &func,
               const Param &param) const;

  /**
   * @brief Return the result of the last @p compute() call.
   * @return Integral value.
   */
  double getSolution() const;

  /**
   * @brief Return the relative error tolerance.
   * @return Tolerance value.
   */
  double getAccuracy() const { return gslIntegrator->getAccuracy(); }

  /**
   * @brief Return the quadrature type.
   * @return Type enum value.
   */
  Type getType() const { return gslIntegrator->getType(); }

private:

  /**
   * @brief Abstract base class for GSL-backed integrators.
   */
  class Base {

  public:

    /**
     * @brief Construct with workspace limit and accuracy.
     * @param type_   Quadrature mode.
     * @param limit_  GSL workspace limit.
     * @param relErr_ Relative error tolerance.
     */
    Base(const Type &type_, const size_t &limit_, const double &relErr_)
        : type(type_),
          limit(limit_),
          relErr(relErr_) {}

    /** @brief Virtual destructor. */
    virtual ~Base() = default;

    /**
     * @brief Return the integral value from the last compute.
     * @return Solution value.
     */
    double getSolution() const { return sol; }

    /** @brief Return the relative error tolerance. */
    double getAccuracy() const { return relErr; }

    /** @brief Return the quadrature type. */
    Type getType() const { return type; }

    /**
     * @brief Compute the integral.
     * @param func  Integrand.
     * @param param Integration parameters.
     */
    virtual void compute(const std::function<double(double)> &func,
                         const Param &param) = 0;

  protected:

    /** @brief GSL function wrapper (populated by subclasses). */
    gsl_function *F;
    /** @brief Quadrature type. */
    const Type type;
    /** @brief GSL workspace limit. */
    const size_t limit;
    /** @brief Relative error tolerance. */
    const double relErr;
    /** @brief Residual absolute error from the last compute. */
    double err;
    /** @brief Integral value from the last compute. */
    double sol;
  };

  /**
   * @brief GSL CQUAD doubly-adaptive quadrature.
   */
  class CQUAD : public Base {
  public:

    /**
     * @brief Construct with a relative error tolerance.
     * @param relErr_ Relative error tolerance.
     */
    explicit CQUAD(const double &relErr_);

    /**
     * @brief Copy constructor.
     * @param other Source CQUAD integrator.
     */
    CQUAD(const CQUAD &other)
        : Integrator1D::CQUAD(other.relErr) {}

    /** @brief Destructor — releases the GSL workspace. */
    ~CQUAD();

    /**
     * @brief Compute the integral.
     * @param func  Integrand.
     * @param param Integration parameters.
     */
    void compute(const std::function<double(double)> &func,
                 const Param &param) override;

  private:

    /** @brief GSL CQUAD workspace. */
    gsl_integration_cquad_workspace *wsp;
    /** @brief Number of integrand evaluations from the last compute. */
    size_t nEvals;
  };

  /**
   * @brief GSL QAWO oscillatory quadrature for Fourier-type integrals.
   */
  class QAWO : public Base {

  public:

    /**
     * @brief Construct with a relative error tolerance.
     * @param relErr_ Relative error tolerance.
     */
    explicit QAWO(const double &relErr_);

    /**
     * @brief Copy constructor.
     * @param other Source QAWO integrator.
     */
    QAWO(const QAWO &other)
        : Integrator1D::QAWO(other.relErr) {}

    /** @brief Destructor — releases the GSL workspaces. */
    ~QAWO();

    /**
     * @brief Compute the integral.
     * @param func  Integrand.
     * @param param Integration parameters.
     */
    void compute(const std::function<double(double)> &func,
                 const Param &param) override;

  private:

    /** @brief Primary GSL integration workspace. */
    gsl_integration_workspace *wsp;
    /** @brief Cycle workspace for the QAWO algorithm. */
    gsl_integration_workspace *wspc;
    /** @brief QAWO table storing the oscillation parameters. */
    gsl_integration_qawo_table *qtab;
  };

  /**
   * @brief GSL QAGS adaptive quadrature with singularity handling.
   */
  class QAGS : public Base {

  public:

    /**
     * @brief Construct with a relative error tolerance.
     * @param relErr_ Relative error tolerance.
     */
    explicit QAGS(const double &relErr_);

    /**
     * @brief Copy constructor.
     * @param other Source QAGS integrator.
     */
    QAGS(const QAGS &other)
        : Integrator1D::QAGS(other.relErr) {}

    /** @brief Destructor — releases the GSL workspace. */
    ~QAGS();

    /**
     * @brief Compute the integral.
     * @param func  Integrand.
     * @param param Integration parameters.
     */
    void compute(const std::function<double(double)> &func,
                 const Param &param) override;

  private:

    /** @brief GSL QAGS workspace. */
    gsl_integration_workspace *wsp;
  };

  /** @brief Pointer to the active GSL integrator backend. */
  std::unique_ptr<Base> gslIntegrator;
};

/**
 * @brief 2D numerical integrator using nested 1D GSL quadrature.
 *
 * Computes @f$\int_{x_{\min}}^{x_{\max}} f_1(x)
 *   \left[ \int_{y_{\min}(x)}^{y_{\max}(x)} f_2(y)\, dy \right] dx@f$
 * by applying an outer 1D integrator over @f$x@f$ and an inner 1D integrator
 * over @f$y@f$, where the inner limits may depend on @f$x@f$.
 */
class Integrator2D {

public:

  /** @brief Alias for the 1D integrator type enum. */
  using Type = Integrator1D::Type;
  /** @brief Alias for the 1D parameter struct. */
  using Param1D = Integrator1D::Param;

  /**
   * @brief 2D integration parameters.
   *
   * Extends the 1D parameter with per-x lower and upper y-limits
   * (which may be constant or x-dependent functions).
   */
  class Param : public Param1D {

  public:

    /** @brief Function type for x-dependent integration limits. */
    using Func = std::function<double(double)>;

    /// @cond INTERNAL
    const Func yMin = [&](const double &x) {
      (void)(x);
      return yMinNum;
    };
    const Func yMax = [&](const double &x) {
      (void)(x);
      return yMaxNum;
    };
    /// @endcond
    /** @brief Optional x-grid for grid-based 2D integration. */
    const std::vector<double> xGrid;

    /**
     * @brief Construct with x-dependent y-limits.
     * @param xMin_ Lower x-limit.
     * @param xMax_ Upper x-limit.
     * @param yMin_ Lower y-limit function of x.
     * @param yMax_ Upper y-limit function of x.
     */
    Param(const double &xMin_,
          const double &xMax_,
          const Func &yMin_,
          const Func &yMax_)
        : Integrator1D::Param(xMin_, xMax_),
          yMin(yMin_),
          yMax(yMax_) {}

    /**
     * @brief Construct with constant y-limits.
     * @param xMin_ Lower x-limit.
     * @param xMax_ Upper x-limit.
     * @param yMin_ Constant lower y-limit.
     * @param yMax_ Constant upper y-limit.
     */
    Param(const double &xMin_,
          const double &xMax_,
          const double &yMin_,
          const double &yMax_)
        : Integrator1D::Param(xMin_, xMax_),
          yMinNum(yMin_),
          yMaxNum(yMax_) {}

    /**
     * @brief Construct for a Fourier-type 2D integral.
     * @param fourierR_ Fourier radius.
     */
    Param(const double &fourierR_)
        : Integrator1D::Param(fourierR_) {}

  private:

    /** @brief Numeric constant for the lower y-limit. */
    const double yMinNum = numUtil::NaN;
    /** @brief Numeric constant for the upper y-limit. */
    const double yMaxNum = numUtil::NaN;
  };

  /**
   * @brief Construct with independent types for the outer and inner
   * integrators.
   * @param type1  Quadrature type for the outer (x) integral.
   * @param type2  Quadrature type for the inner (y) integral.
   * @param relErr Relative error tolerance for both integrators.
   */
  Integrator2D(const Type &type1, const Type &type2, const double &relErr)
      : itg1(type1, relErr),
        itg2(type2, relErr) {}

  /**
   * @brief Construct with the same quadrature type for both levels.
   * @param type   Quadrature type.
   * @param relErr Relative error tolerance.
   */
  Integrator2D(const Type &type, const double &relErr)
      : Integrator2D(type, type, relErr) {}

  /**
   * @brief Construct using the default (CQUAD) quadrature type.
   * @param relErr Relative error tolerance.
   */
  Integrator2D(const double &relErr)
      : Integrator2D(Type::DEFAULT, relErr) {}

  /**
   * @brief Compute the double integral.
   * @param func1  Outer integrand @f$f_1(x)@f$.
   * @param func2  Inner integrand @f$f_2(y)@f$.
   * @param param  2D integration parameters.
   * @param xGrid  Grid for grid-based outer integration.
   */
  void compute(const std::function<double(double)> &func1,
               const std::function<double(double)> &func2,
               const Param &param,
               const std::vector<double> &xGrid);

  /**
   * @brief Return the current outer-integration variable.
   * @return Most recent value of @f$x@f$ used in the inner integral.
   */
  double getX() const { return x; };

  /**
   * @brief Return the result of the last @p compute() call.
   * @return Double-integral value.
   */
  double getSolution() const { return sol; };

private:

  /** @brief Outer (level 1) integrator over @f$x@f$. */
  Integrator1D itg1;
  /** @brief Inner (level 2) integrator over @f$y@f$. */
  Integrator1D itg2;
  /** @brief Temporary storage for the current outer variable @f$x@f$. */
  double x;
  /** @brief Integral result from the last @p compute() call. */
  double sol;
};

#endif
