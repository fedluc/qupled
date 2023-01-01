#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <vector>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

using namespace std;

// Wrapper for gsl_function objects
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

// Root solver
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


// Integrator for 1D integrals
class Integrator1D {

private:

  // Function to integrate
  gsl_function *F;
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
  = gsl_integration_cquad_workspace_alloc(100);
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
  Integrator1D(const double relErr_) : relErr(relErr_) {;};
  // Destructor
  ~Integrator1D(){
    gsl_integration_cquad_workspace_free(wsp);
  }
  void compute(const function<double(double)> func,
	       const double xMin,
	       const double xMax);
  double getSolution() { return sol; };
  
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
  void compute(const function<double(double)> func);
  double getSolution() { return sol; };
  
};

// Data interpolator
class Interpolator {

private:

  // Data to interpolate
  const vector<double> xData;
  const vector<double> yData;
  // Spline
  gsl_spline *spline;
  // Accelerator
  gsl_interp_accel *acc;
  
  
public:

  // Constructor
  Interpolator(vector<double> xData_, vector<double> yData_)
    : xData(xData_), yData(yData_) {
    assert(xData.size() == yData.size());
    spline = gsl_spline_alloc(gsl_interp_cspline, xData.size());
    acc = gsl_interp_accel_alloc();
    gsl_spline_init(spline, &xData[0], &yData[0], xData.size());
  };
  // Destructor
  ~Interpolator(){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  // Evaluate
  double eval(double x) const {
    return gsl_spline_eval(spline, x, acc);
  };
  
};



#endif
