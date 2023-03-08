#include "numerics.hpp"

// -----------------------------------------------------------------
// Interpolator class
// -----------------------------------------------------------------

// Constructors
Interpolator::Interpolator(const double &x,
			   const double &y,
			   const size_t n_) {
  n = n_;
  spline = gsl_spline_alloc(gsl_interp_cspline, n);
  acc = gsl_interp_accel_alloc();
  gsl_spline_init(spline, &x, &y, n);
};

Interpolator::Interpolator(const Interpolator &it) {
  if (spline) gsl_spline_free(spline);
  if (acc) gsl_interp_accel_free(acc);
  n = it.n;
  spline = gsl_spline_alloc(gsl_interp_cspline, n);
  acc = gsl_interp_accel_alloc();
  *spline = *it.spline;
  *acc = *it.acc;
};

Interpolator::Interpolator() {
  n = 0;
  spline = gsl_spline_alloc(gsl_interp_cspline, n);
  acc = gsl_interp_accel_alloc();
};

// Destructor
Interpolator::~Interpolator(){
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

// Evaluate interpolation
double Interpolator::eval(const double x) const {
  return gsl_spline_eval(spline, x, acc);
};

// -----------------------------------------------------------------
// Interpolator2D class
// -----------------------------------------------------------------

// Constructors
Interpolator2D::Interpolator2D() {
  spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, 1, 1);
  xacc = gsl_interp_accel_alloc();
  yacc = gsl_interp_accel_alloc();
};

    
Interpolator2D::Interpolator2D(const double &x,
			       const double &y,
			       const double &z,
			       const int nx_,
			       const int ny_) {
  nx = nx_;
  ny = ny_;
  spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, nx, ny);
  xacc = gsl_interp_accel_alloc();
  yacc = gsl_interp_accel_alloc();
  // Ensure that z is stored in the correct order
  double *za = (double*)malloc(nx * ny * sizeof(double));
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      gsl_spline2d_set(spline, za, i, j, *(&z + j + i*ny)); 
    }
  }
  gsl_spline2d_init(spline, &x, &y, za, nx, ny);
  free(za);
};
  
Interpolator2D::Interpolator2D(const Interpolator2D &it) {
  if (spline) gsl_spline2d_free(spline);
  if (xacc) gsl_interp_accel_free(xacc);
  if (yacc) gsl_interp_accel_free(yacc);
  nx = it.nx;
  ny = it.ny;
  spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, nx, ny);
  xacc = gsl_interp_accel_alloc();
  yacc = gsl_interp_accel_alloc();
  *spline = *it.spline;
  *xacc = *it.xacc;
  *yacc = *it.yacc;
};

// Destructor
Interpolator2D::~Interpolator2D(){
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
}

// Evaluate
double Interpolator2D::eval(const double x, const double y) const {
  return gsl_spline2d_eval(spline, x, y, xacc, yacc);
};

// -----------------------------------------------------------------
// RootSolver class
// -----------------------------------------------------------------

void RootSolver::solve(const function<double(double)> func,
		       const vector<double> guess){
  // Set up function
  GslFunctionWrap<decltype(func)> Fp(func);
  F = static_cast<gsl_function*>(&Fp);
  // Set up solver
  gsl_root_fsolver_set(rs, F, guess.at(0), guess.at(1));
  // Call solver
  do{
    status = gsl_root_fsolver_iterate(rs);
    sol = gsl_root_fsolver_root(rs);
    double solLo = gsl_root_fsolver_x_lower(rs);
    double solHi = gsl_root_fsolver_x_upper(rs);
    status = gsl_root_test_interval(solLo, solHi, 0, relErr);
    iter++;
  } while (status == GSL_CONTINUE && iter < maxIter);
}

// -----------------------------------------------------------------
// Integrator1D class
// -----------------------------------------------------------------

void Integrator1D::compute(const function<double(double)> func,
			   const double xMin,
			   const double xMax){
  // Set up function
  GslFunctionWrap<decltype(func)> Fp(func);
  F = static_cast<gsl_function*>(&Fp);
  // Integrate
  gsl_integration_cquad(F, xMin, xMax, 
			0.0, relErr, 
			wsp, &sol,
			&err, &nEvals);
}

// -----------------------------------------------------------------
// Integrator1DFourier class
// -----------------------------------------------------------------

void Integrator1DFourier::compute(const function<double(double)> func){
  // Set up function
  GslFunctionWrap<decltype(func)> Fp(func);
  F = static_cast<gsl_function*>(&Fp);
  // Set wave-vector
  gsl_integration_qawo_table_set(qtab, r, 1.0, GSL_INTEG_SINE);
  // Integrate
  gsl_integration_qawf(F, 0.0, relErr, 
		       limit, wsp, wspc,
		       qtab, &sol, &err);
}

// -----------------------------------------------------------------
// Integrator2D class
// -----------------------------------------------------------------

void Integrator2D::compute(const function<double(double)> func1,
			   const function<double(double)> func2,
			   const double xMin,
			   const double xMax,
			   const function<double(double)> yMin,
			   const function<double(double)> yMax){
  // Level 2 integration
  auto func = [&](double x_)->double {
    x = x_;
    itg2.compute(func2, yMin(x_), yMax(x_));
    return func1(x_) * itg2.getSolution();
  };
  itg1.compute(func, xMin, xMax);
  sol = itg1.getSolution();
}

