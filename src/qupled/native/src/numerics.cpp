#include "numerics.hpp"
#include "mpi_util.hpp"
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std;
using namespace GslWrappers;
using namespace MPIUtil;

// -----------------------------------------------------------------
// C++ wrappers to GSL objects
// -----------------------------------------------------------------

template <typename Func, typename... Args>
void GslWrappers::callGSLFunction(Func &&gslFunction, Args &&...args) {
  int status = gslFunction(std::forward<Args>(args)...);
  if (status) {
    throwError("GSL error: " + std::to_string(status) + ", "
               + std::string(gsl_strerror(status)));
  }
}

template <typename Ptr, typename Func, typename... Args>
void GslWrappers::callGSLAlloc(Ptr &ptr, Func &&gslFunction, Args &&...args) {
  ptr = gslFunction(std::forward<Args>(args)...);
  if (!ptr) { throwError("GSL error: allocation error"); }
}

// -----------------------------------------------------------------
// Wrappers to GSL special functions
// -----------------------------------------------------------------

double SpecialFunctions::fermiDirac(const double &x) {
  gsl_sf_result gamma, fd;
  callGSLFunction(gsl_sf_gamma_e, 1.5, &gamma);
  callGSLFunction(gsl_sf_fermi_dirac_half_e, x, &fd);
  return gamma.val * fd.val;
}

// -----------------------------------------------------------------
// Interpolator class
// -----------------------------------------------------------------

// Constructors
Interpolator1D::Interpolator1D() {
  n = 0;
  spline = nullptr;
  acc = nullptr;
}

Interpolator1D::Interpolator1D(const vector<double> &x, const vector<double> &y)
    : Interpolator1D() {
  assert(x.size() == y.size());
  setup(x[0], y[0], x.size());
}

Interpolator1D::Interpolator1D(const double &x,
                               const double &y,
                               const size_t n_)
    : Interpolator1D() {
  setup(x, y, n_);
}

// Destructor
Interpolator1D::~Interpolator1D() {
  if (spline) gsl_spline_free(spline);
  if (acc) gsl_interp_accel_free(acc);
}

// Setup interpolator
void Interpolator1D::setup(const double &x, const double &y, const size_t n_) {
  n = n_;
  if (!isValid()) {
    n = 0;
    return;
  }
  cutoff = *(&x + n - 1);
  callGSLAlloc(spline, gsl_spline_alloc, TYPE, n);
  callGSLAlloc(acc, gsl_interp_accel_alloc);
  callGSLFunction(gsl_spline_init, spline, &x, &y, n);
}

// Check if the interpolator can be setup correctly
bool Interpolator1D::isValid() const {
  return n >= gsl_interp_type_min_size(TYPE);
}

// Reset existing interpolator
void Interpolator1D::reset(const double &x, const double &y, const size_t n_) {
  if (spline) gsl_spline_free(spline);
  if (acc) gsl_interp_accel_free(acc);
  setup(x, y, n_);
}

// Evaluate interpolation
double Interpolator1D::eval(const double &x) const {
  double out;
  callGSLFunction(
      gsl_spline_eval_e, spline, (x < cutoff) ? x : cutoff, acc, &out);
  return out;
}

// -----------------------------------------------------------------
// Interpolator2D class
// -----------------------------------------------------------------

// Constructors
Interpolator2D::Interpolator2D() {
  nx = 0;
  ny = 0;
  spline = nullptr;
  xacc = nullptr;
  yacc = nullptr;
}

Interpolator2D::Interpolator2D(const double &x,
                               const double &y,
                               const double &z,
                               const int nx_,
                               const int ny_)
    : Interpolator2D() {
  setup(x, y, z, nx_, ny_);
}

// Destructor
Interpolator2D::~Interpolator2D() {
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
}

// Check if the interpolator can be setup correctly
bool Interpolator2D::isValid() const {
  return nx + ny >= gsl_interp2d_type_min_size(TYPE);
}

// Setup interpolator
void Interpolator2D::setup(const double &x,
                           const double &y,
                           const double &z,
                           const int nx_,
                           const int ny_) {
  nx = nx_;
  ny = ny_;
  if (!isValid()) {
    nx = 0;
    ny = 0;
    return;
  }
  callGSLAlloc(spline, gsl_spline2d_alloc, TYPE, nx, ny);
  callGSLAlloc(xacc, gsl_interp_accel_alloc);
  callGSLAlloc(yacc, gsl_interp_accel_alloc);
  // Ensure that z is stored in the correct order
  double *za = (double *)malloc(nx * ny * sizeof(double));
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      callGSLFunction(gsl_spline2d_set, spline, za, i, j, *(&z + j + i * ny));
    }
  }
  callGSLFunction(gsl_spline2d_init, spline, &x, &y, za, nx, ny);
  free(za);
}

// Reset existing interpolator
void Interpolator2D::reset(const double &x,
                           const double &y,
                           const double &z,
                           const int nx_,
                           const int ny_) {
  if (spline) gsl_spline2d_free(spline);
  if (xacc) gsl_interp_accel_free(xacc);
  if (yacc) gsl_interp_accel_free(yacc);
  setup(x, y, z, nx_, ny_);
}

// Evaluate interpolation
double Interpolator2D::eval(const double &x, const double &y) const {
  double out;
  callGSLFunction(gsl_spline2d_eval_e, spline, x, y, xacc, yacc, &out);
  return out;
}

// -----------------------------------------------------------------
// BrentRootSolver class
// -----------------------------------------------------------------

// Constructor
BrentRootSolver::BrentRootSolver()
    : rst(gsl_root_fsolver_brent) {
  callGSLAlloc(rs, gsl_root_fsolver_alloc, rst);
}

// Destructor
BrentRootSolver::~BrentRootSolver() { gsl_root_fsolver_free(rs); }

// Invoke root solver
void BrentRootSolver::solve(const function<double(double)> &func,
                            const vector<double> &guess) {
  // Set up function
  GslFunctionWrap<decltype(func)> Fp(func);
  F = static_cast<gsl_function *>(&Fp);
  // Set up solver
  callGSLFunction(gsl_root_fsolver_set, rs, F, guess.at(0), guess.at(1));
  // Call solver
  do {
    callGSLFunction(gsl_root_fsolver_iterate, rs);
    sol = gsl_root_fsolver_root(rs);
    double solLo = gsl_root_fsolver_x_lower(rs);
    double solHi = gsl_root_fsolver_x_upper(rs);
    status = gsl_root_test_interval(solLo, solHi, 0, relErr);
    iter++;
  } while (status == GSL_CONTINUE && iter < maxIter);
  // Check if the solver managed to find a solution
  if (status != GSL_SUCCESS) {
    throwError("The brent root solver "
               "did not converge to the desired accuracy.");
  }
}

// -----------------------------------------------------------------
// SecantSolver class
// -----------------------------------------------------------------

void SecantSolver::solve(const function<double(double)> &func,
                         const vector<double> &guess) {
  // Set up solver
  double x0 = guess.at(0);
  double x1 = guess.at(1);
  double fx0;
  double fx1 = func(x0);
  // Call solver
  do {
    fx0 = fx1;
    fx1 = func(x1);
    sol = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
    if (abs(sol - x1) < abs(sol) * relErr) { status = GSL_SUCCESS; }
    x0 = x1;
    x1 = sol;
    iter++;
  } while (status == GSL_CONTINUE && iter < maxIter);
  // Check if the solver managed to find a solution
  if (status != GSL_SUCCESS) {
    throwError("The secant root solver "
               "did not converge to the desired accuracy.");
  }
}

// -----------------------------------------------------------------
// Integrator1D class
// -----------------------------------------------------------------

// Constructor
Integrator1D::Integrator1D(const Type &type, const double &relErr) {
  switch (type) {
  case Type::DEFAULT: gslIntegrator = make_unique<CQUAD>(relErr); break;
  case Type::FOURIER: gslIntegrator = make_unique<QAWO>(relErr); break;
  case Type::SINGULAR: gslIntegrator = make_unique<QAGS>(relErr); break;
  default: throwError("Invalid integrator type");
  }
}

// Compute integral
void Integrator1D::compute(const std::function<double(double)> &func,
                           const Param &param) const {
  gslIntegrator->compute(func, param);
}

// Getters
double Integrator1D::getSolution() const {
  return gslIntegrator->getSolution();
}

// -----------------------------------------------------------------
// Integrator1D::CQUAD class
// -----------------------------------------------------------------

// Constructor
Integrator1D::CQUAD::CQUAD(const double &relErr_)
    : Integrator1D::Base(Type::DEFAULT, 100, relErr_) {
  callGSLAlloc(wsp, gsl_integration_cquad_workspace_alloc, limit);
}

// Destructor
Integrator1D::CQUAD::~CQUAD() { gsl_integration_cquad_workspace_free(wsp); }

// Compute integral
void Integrator1D::CQUAD::compute(const function<double(double)> &func,
                                  const Param &param) {
  // Check parameter validity
  if (isnan(param.xMin) || isnan(param.xMax)) {
    throwError("Integration limits were not set correctly");
  }
  // Set up function
  GslFunctionWrap<decltype(func)> Fp(func);
  F = static_cast<gsl_function *>(&Fp);
  // Integrate
  callGSLFunction(gsl_integration_cquad,
                  F,
                  param.xMin,
                  param.xMax,
                  0.0,
                  relErr,
                  wsp,
                  &sol,
                  &err,
                  &nEvals);
}

// -----------------------------------------------------------------
// Integrator1D::QAWO class
// -----------------------------------------------------------------

// Constructor
Integrator1D::QAWO::QAWO(const double &relErr_)
    : Integrator1D::Base(Type::FOURIER, 1000, relErr_) {
  callGSLAlloc(wsp, gsl_integration_workspace_alloc, limit);
  callGSLAlloc(wspc, gsl_integration_workspace_alloc, limit);
  callGSLAlloc(
      qtab, gsl_integration_qawo_table_alloc, 0.0, 1.0, GSL_INTEG_SINE, limit);
}

// Destructor
Integrator1D::QAWO::~QAWO() {
  gsl_integration_workspace_free(wsp);
  gsl_integration_workspace_free(wspc);
  gsl_integration_qawo_table_free(qtab);
}

// Compute integral
void Integrator1D::QAWO::compute(const function<double(double)> &func,
                                 const Param &param) {
  // Check parameter validity
  if (isnan(param.fourierR)) {
    throwError("Integration parameters were not set correctly");
  }
  // Set up function
  GslFunctionWrap<decltype(func)> Fp(func);
  F = static_cast<gsl_function *>(&Fp);
  // Set wave-vector
  callGSLFunction(gsl_integration_qawo_table_set,
                  qtab,
                  param.fourierR,
                  1.0,
                  GSL_INTEG_SINE);
  // Integrate
  callGSLFunction(
      gsl_integration_qawf, F, 0.0, relErr, limit, wsp, wspc, qtab, &sol, &err);
}

// -----------------------------------------------------------------
// Integrator1D::QAGS class
// -----------------------------------------------------------------

// Constructor
Integrator1D::QAGS::QAGS(const double &relErr_)
    : Integrator1D::Base(Type::SINGULAR, 1000, relErr_) {
  callGSLAlloc(wsp, gsl_integration_workspace_alloc, limit);
}

// Destructor
Integrator1D::QAGS::~QAGS() { gsl_integration_workspace_free(wsp); }

// Compute integral
void Integrator1D::QAGS::compute(const function<double(double)> &func,
                                 const Param &param) {
  // Check parameter validity
  if (isnan(param.xMin) || isnan(param.xMax)) {
    throwError("Integration limits were not set correctly");
  }
  // Set up function
  GslFunctionWrap<decltype(func)> Fp(func);
  F = static_cast<gsl_function *>(&Fp);
  // Integrate
  callGSLFunction(gsl_integration_qags,
                  F,
                  param.xMin,
                  param.xMax,
                  0.0,
                  relErr,
                  limit,
                  wsp,
                  &sol,
                  &err);
}

// -----------------------------------------------------------------
// Integrator2D class
// -----------------------------------------------------------------

// Compute integral
void Integrator2D::compute(const function<double(double)> &func1,
                           const function<double(double)> &func2,
                           const Param &param,
                           const vector<double> &xGrid) {
  const int nx = xGrid.size();
  function<double(double)> func;
  Interpolator1D itp;
  auto getParam1D = [&](const double &x) {
    const bool isFourier = !isnan(param.fourierR);
    if (isFourier) { return Param1D(param.fourierR); }
    return Param1D(param.yMin(x), param.yMax(x));
  };
  if (nx > 0) {
    // Level 2 integration (only evaluated at the points in xGrid)
    vector<double> sol2(nx);
    for (int i = 0; i < nx; ++i) {
      x = xGrid[i];
      itg2.compute(func2, getParam1D(x));
      sol2[i] = itg2.getSolution();
    }
    itp.reset(xGrid[0], sol2[0], nx);
    func = [&](const double &x_) -> double { return func1(x_) * itp.eval(x_); };
  } else {
    // Level 2 integration (evaluated at arbitrary points)
    func = [&](const double &x_) -> double {
      x = x_;
      itg2.compute(func2, getParam1D(x));
      return func1(x_) * itg2.getSolution();
    };
  }
  // Level 1 integration
  itg1.compute(func, param);
  sol = itg1.getSolution();
}
