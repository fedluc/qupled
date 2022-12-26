#include <gsl/gsl_errno.h>
#include "numerics.hpp"

// Run root solver
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

// Compute integrals 
void Integrator1D::compute(const function<double(double)> func,
			   const double xMin, double xMax){
  // Set up function
  GslFunctionWrap<decltype(func)> Fp(func);
  F = static_cast<gsl_function*>(&Fp);
  // Integrate
  gsl_integration_cquad(F, xMin, xMax, 
			0.0, relErr, 
			wsp, &sol,
			&err, &nEvals);
}

// Element-wise sum between two vectors
vector<double> vecUtil::sum(const vector<double> &v1,
			const vector<double> &v2) {
  assert(v1.size() == v2.size());
  vector<double> res;
  transform(v1.begin(), v1.end(),
	    v2.begin(), back_inserter(res),
	    plus<double>());
  return res;
}

// Element-wise difference between two vectors
vector<double> vecUtil::diff(const vector<double> &v1,
			     const vector<double> &v2) {
  assert(v1.size() == v2.size());
  vector<double> res;
  transform(v1.begin(), v1.end(),
	    v2.begin(), back_inserter(res),
	    minus<double>());
  return res;
}

// Root square difference between two vectors
double vecUtil::rms(const vector<double> &v1,
		    const vector<double> &v2,
		    const bool normalize) {
  const vector<double> tmp = diff(v1,v2);
  double rms = inner_product(tmp.begin(), tmp.end(), tmp.begin(), 0.0);
  if (normalize) rms /= tmp.size();
  return sqrt(rms);
}

// Element-wise multiplication of a vector and a scalar
vector<double> vecUtil::mult(const vector<double> &v,
			     const double a) {
  vector<double> res = v;
  transform(res.begin(), res.end(), res.begin(), [&a](double c){return c*a;});
  return res;
}
