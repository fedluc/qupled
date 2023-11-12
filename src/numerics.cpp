#include "numerics.hpp"

// -----------------------------------------------------------------
// Interpolator class
// -----------------------------------------------------------------

// Constructors
Interpolator1D::Interpolator1D(const vector<double> &x,
			       const vector<double> &y) {
  assert(x.size() == y.size());
  setup(x[0], y[0], x.size());
}

Interpolator1D::Interpolator1D(const double &x,
			       const double &y,
			       const size_t n_) {
  setup(x, y, n_);
}

Interpolator1D::Interpolator1D() {
  n = 0;
  spline = nullptr;
  acc = nullptr;
}

// Destructor
Interpolator1D::~Interpolator1D(){
  if (spline) gsl_spline_free(spline);
  if (acc) gsl_interp_accel_free(acc);
}

// Setup interpolator
void Interpolator1D::setup(const double &x,
			   const double &y,
			   const size_t n_) {
  n = n_;
  cutoff = *(&x + n - 1);
  spline = gsl_spline_alloc(gsl_interp_cspline, n);
  acc = gsl_interp_accel_alloc();
  gsl_spline_init(spline, &x, &y, n);
}

// Reset existing interpolator
void Interpolator1D::reset(const double &x,
			   const double &y,
			   const size_t n_) {
  if (spline) gsl_spline_free(spline);
  if (acc) gsl_interp_accel_free(acc);
  setup(x, y, n_);
}

// Evaluate interpolation
double Interpolator1D::eval(const double x) const {
  return gsl_spline_eval(spline, (x < cutoff) ? x : cutoff, acc);
}

// -----------------------------------------------------------------
// Interpolator2D class
// -----------------------------------------------------------------

// Constructors    
Interpolator2D::Interpolator2D(const double &x,
			       const double &y,
			       const double &z,
			       const int nx_,
			       const int ny_) {
  setup(x, y, z, nx_, ny_);
}

Interpolator2D::Interpolator2D() {
  nx = 0;
  ny = 0;
  spline = nullptr;
  xacc = nullptr;
  yacc = nullptr;
}

// Destructor
Interpolator2D::~Interpolator2D(){
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
}

// Setup interpolator
void Interpolator2D::setup(const double &x,
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
double Interpolator2D::eval(const double x, const double y) const {
  return gsl_spline2d_eval(spline, x, y, xacc, yacc);
};

// -----------------------------------------------------------------
// BrentRootSolver class
// -----------------------------------------------------------------

void BrentRootSolver::solve(const function<double(double)> func,
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

void SecantSolver::solve(const function<double(double)> func,
			 const vector<double> guess){
  // Set up solver
  double x0 = guess.at(0);
  double x1 = guess.at(1);
  double fx0;
  double fx1 = func(x0);
  status = GSL_CONTINUE;
  // Call solver
  do{
    fx0 = fx1;
    fx1 = func(x1);
    sol = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
    if ( abs(sol - x1) < abs(sol) * relErr) { status = GSL_SUCCESS; }
    x0 = x1;
    x1 = sol;
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
			   const function<double(double)> yMax,
			   const vector<double> &xGrid){
  const int nx = xGrid.size();
  function<double(double)> func;
  Interpolator1D itp;
  if (nx > 0) {
    // Level 2 integration (only evaluated at the points in xGrid)
    vector<double> sol2(nx);
    for (int i = 0; i < nx; ++i) {
      x = xGrid[i];
      itg2.compute(func2, yMin(x), yMax(x));
      sol2[i] = itg2.getSolution();
    }
    itp.reset(xGrid[0], sol2[0], nx);
    func = [&](double x_)->double {
      return func1(x_) * itp.eval(x_);
    };
  }
  else {
    // Level 2 integration (evaluated at arbitrary points) 
    func = [&](double x_)->double {
      x = x_;
      itg2.compute(func2, yMin(x_), yMax(x_));
      return func1(x_) * itg2.getSolution();
    };
  }
  // Level 1 integration
  itg1.compute(func, xMin, xMax);
  sol = itg1.getSolution();
}


void Integrator2D::compute(const function<double(double)> func1,
			   const function<double(double)> func2,
			   const double xMin,
			   const double xMax,
			   const function<double()> yMin,
			   const function<double()> yMax,
			   const vector<double> &xGrid){
  // Wrappers around yMin and yMax to avoid compiler warnings
  auto yMinTmp = [&](double x){(void)(x); return yMin(); };
  auto yMaxTmp = [&](double x){(void)(x); return yMax(); };
  compute(func1, func2, xMin, xMax, yMinTmp, yMaxTmp, xGrid);
}
