#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <vector>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

using namespace std;

// -----------------------------------------------------------------
// Classes to interpolate data
// -----------------------------------------------------------------

// Interpolator for 1D data
class Interpolator1D {

private:

  // Spline
  gsl_spline *spline;
  // Accelerator
  gsl_interp_accel *acc;
  // Size
  size_t n;
  // 
  // Setup interpolator
  void setup(const double &x,
	     const double &y,
	     const size_t sz_);
  
public:

  // Constructor
  Interpolator1D(const vector<double> &x,
		 const vector<double> &y);
  Interpolator1D(const double &x,
		 const double &y,
		 const size_t sz_);
  Interpolator1D();
  // Destructor
  ~Interpolator1D();
  // Reset
  void reset(const double &x,
	     const double &y,
	     const size_t sz_);
  // Evaluate
  double eval(const double x) const;
  
};

// Interpolator for 2D data
class Interpolator2D {

private:

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
	     const int szx_,
	     const int szy_);
public:

  // Constructor
  Interpolator2D(const double &x,
		 const double &y,
		 const double &z,
		 const int szx_,
		 const int szy_);
  Interpolator2D(const Interpolator2D &it);
  Interpolator2D();
  // Destructor
  ~Interpolator2D();
  // Reset
  void reset(const double &x,
	     const double &y,
	     const double &z,
	     const int szx_,
	     const int szy_);
  // Evaluate
  double eval(const double x, const double y) const;
};

// -----------------------------------------------------------------
// C++ wrappers to gsl_function objects
// -----------------------------------------------------------------

template< typename T >
class GslFunctionWrap : public gsl_function {

private:

  const T& func;
  static double invoke(double x, void *params) {
    return static_cast<GslFunctionWrap*>(params)->func(x);
  }
  
public:
  
  GslFunctionWrap(const T& func_) : func(func_) {
    function = &GslFunctionWrap::invoke;
    params=this;
  }
  
};

  
// -----------------------------------------------------------------
// Classes to find roots of equations
// -----------------------------------------------------------------

class RootSolver {

private:

  // Function to solve
  gsl_function *F;
  // Type of solver
  const gsl_root_fsolver_type *rst = gsl_root_fsolver_bisection;
  // Solver
  gsl_root_fsolver *rs = gsl_root_fsolver_alloc(rst);
  // Accuracy
  const double relErr = 1e-10;
  // Iterations
  const int maxIter = 1000;
  int iter = 0;
  // Solver status
  int status;
  // Solution
  double sol;
  
public:

  RootSolver() {;};
  RootSolver(const double relErr_) : relErr(relErr_) {;};
  RootSolver(const int maxIter_) : maxIter(maxIter_) {;};
  RootSolver(const double relErr_, const int maxIter_) :
    relErr(relErr_), maxIter(maxIter_) {;};
  void solve(const function<double(double)> func,
	     const vector<double> guess);
  bool success() const { return status == GSL_SUCCESS; };
  double getSolution() const { return sol; };
  
};

// -----------------------------------------------------------------
// Classes to compute integrals
// -----------------------------------------------------------------

// Integrator for 1D integrals
class Integrator1D {

private:

  // Function to integrate
  gsl_function *F;
  // Integration workspace limit
  const size_t limit = 100;
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
  = gsl_integration_cquad_workspace_alloc(limit);
  // Accuracy
  const double relErr = 1e-5;
  // Residual error
  double err;
  // Number of evaluations
  size_t nEvals;
  // Solution
  double sol;
  
public:

  // Constructors
  Integrator1D() {;};
  // Destructor
  ~Integrator1D(){
    gsl_integration_cquad_workspace_free(wsp);
  }
  // Compute integral
  void compute(const function<double(double)> func,
	       const double xMin,
	       const double xMax);
  // Getters
  double getSolution() const { return sol; };
  
};

// Integrator for 2D integrals
class Integrator2D {

private:

  // Level 1 integrator (outermost integral)
  Integrator1D itg1;
  // Level 2 integrator
  Integrator1D itg2;
  // Temporary variable for level 2 integration
  double x;
  // Solution
  double sol;
  
public:

  // Constructors
  Integrator2D() {;};
  // Compute integral
  void compute(const function<double(double)> func1,
	       const function<double(double)> func2,
	       const double xMin,
	       const double xMax,
	       const function<double(double)> yMin,
	       const function<double(double)> yMax);
  // Getters
  double getX() const { return x; };
  double getSolution() const { return sol; };
};


// Integrator for 1D integrals of Fourier type
class Integrator1DFourier {

private:

  // Function to integrate
  gsl_function *F;
  // Integration workspace limit
  const size_t limit = 1000;
  // Integration workspace
  gsl_integration_workspace *wsp
  = gsl_integration_workspace_alloc(limit);
  gsl_integration_workspace *wspc
  = gsl_integration_workspace_alloc(limit);
  gsl_integration_qawo_table *qtab
  = gsl_integration_qawo_table_alloc(0.0, 1.0, GSL_INTEG_SINE, limit);
  // Spatial position
  double r;
  const double relErr = 1e-5;
  // Residual error
  double err;
  // Solution
  double sol;
  
public:

  // Constructors
  Integrator1DFourier(const double r_) : r(r_) {;};
  Integrator1DFourier(const double r_, const double relErr_)
    : r(r_), relErr(relErr_) {;};
  // Set spatial position (to re-use the integrator for different r)
  void setR(const double r_) {r = r_;};
  // Destructor
  ~Integrator1DFourier(){
    gsl_integration_workspace_free(wsp);
    gsl_integration_workspace_free(wspc); 
    gsl_integration_qawo_table_free(qtab);
  }
  // Compute integral
  void compute(const function<double(double)> func);
  // Getters
  double getSolution() const { return sol; };
  
};

#endif
