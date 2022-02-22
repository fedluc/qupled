#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "chemical_potential.h"

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

struct nc_params {

  double Theta;

};

double normalization_condition(double mu, void *pp) {

  struct nc_params *params = (struct nc_params*)pp;
  double Theta = (params->Theta);

  return gsl_sf_gamma(1.5)*gsl_sf_fermi_dirac_half(mu) 
    - 2.0/(3.0*pow(Theta, 3.0/2.0));

}

double compute_chemical_potential(input in) {
  
  // Variables
  double mu_lo = in.mu_lo;
  double mu_hi = in.mu_hi;
  double mu;
  int status, iter;

  // Set-up function
  gsl_function ff_root;
  ff_root.function = &normalization_condition;
  struct nc_params ncp = {in.Theta};
  ff_root.params = &ncp;

  // Set-up root-solver
  const gsl_root_fsolver_type * rst = gsl_root_fsolver_bisection;
  gsl_root_fsolver * rs = gsl_root_fsolver_alloc(rst);
  gsl_root_fsolver_set(rs, &ff_root, mu_lo, mu_hi);

  // Solve normalization condition to find chemical potential
  iter = 0;
  do
  {
    
    // Solver iteration
    status = gsl_root_fsolver_iterate(rs);

    // Get solver status
    mu = gsl_root_fsolver_root(rs);
    mu_lo = gsl_root_fsolver_x_lower(rs);
    mu_hi = gsl_root_fsolver_x_upper(rs);
    status = gsl_root_test_interval (mu_lo, mu_hi,
				     0, ROOTMIN_REL_ERR);
    
    // Update iteration counter
    iter++;

  }
  while (status == GSL_CONTINUE && iter < ROOTMIN_MAX_ITER);

  // Free memory
  gsl_root_fsolver_free(rs);

  // Output
  return mu;

}

