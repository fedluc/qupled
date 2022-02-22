#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "solvers.h"
#include "chemical_potential.h"
#include "stls.h"


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS EQUATIONS
// -------------------------------------------------------------------

void solve_stls(input in, bool verbose) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *GG_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;

  // Allocate arrays
  alloc_stls_arrays(in, &xx, &phi, &GG, &GG_new, &SS, &SSHF);

  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_stls_arrays(&in, xx, phi, SSHF, verbose);
  
  // Initial guess
  if (strcmp(in.stls_guess_file,"NO_FILE")==0){
    for (int ii=0; ii < in.nx; ii++) {
      GG[ii] = 0.0; // Static local field correction
      GG_new[ii] = 1.0;
    }
    compute_ssf_stls(SS, SSHF, GG, phi, xx, in); // Static structure factor
  }
  else {
    read_guess_stls(SS, GG, in); // Read from file
  }


  // Iterative procedure
  /* if (verbose) printf("SSF and SLFC calculation...\n"); */
  /* double iter_err = 1.0; */
  /* int iter_counter = 0; */
  /* while (iter_counter < in.nIter && iter_err > in.err_min_iter ) { */
    
  /*   // Start timing */
  /*   double tic = omp_get_wtime(); */
    
  /*   // Update SSF */
  /*   compute_ssf_stls(SS, SSHF, GG, phi, xx, in); */

  /*   // Update SLFC */
  /*   compute_slfc(GG_new, SS, xx, in); */

  /*   // Update diagnostic */
  /*   iter_err = 0.0; */
  /*   iter_counter++; */
  /*   for (int ii=0; ii<in.nx; ii++) { */
  /*     iter_err += (GG_new[ii] - GG[ii]) * (GG_new[ii] - GG[ii]); */
  /*     GG[ii] = in.a_mix*GG_new[ii] + (1-in.a_mix)*GG[ii]; */
  /*   } */
  /*   iter_err = sqrt(iter_err);   */
   
  /*   // End timing */
  /*   double toc = omp_get_wtime(); */
    
  /*   // Print diagnostic */
  /*   if (verbose) { */
  /*     printf("--- iteration %d ---\n", iter_counter); */
  /*     printf("Elapsed time: %f seconds\n", toc - tic); */
  /*     printf("Residual error: %.5e\n", iter_err); */
  /*     fflush(stdout); */
  /*   } */
  /* } */
  if (verbose) printf("Done.\n");
  
  // Internal energy
  if (verbose) printf("Internal energy: %.10f\n",compute_internal_energy(SS, xx, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text_stls(SS, GG, phi, SSHF, xx, in);
  write_guess_stls(SS, GG, in); 
  if (verbose) printf("Done.\n");

  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);

 
}

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays(input in, double **xx, double **phi, 
		       double **GG, double **GG_new, 
		       double **SS, double **SSHF){

  *xx = malloc( sizeof(double) * in.nx);
  *phi = malloc( sizeof(double) * in.nx * in.nl);
  *SSHF = malloc( sizeof(double) * in.nx);
  *GG = malloc( sizeof(double) * in.nx);
  *GG_new = malloc( sizeof(double) * in.nx);
  *SS = malloc( sizeof(double) * in.nx);
  
}

void free_stls_arrays(double *xx, double *phi, double *GG, 
		      double *GG_new, double *SS,
		      double *SSHF){

  free(xx);
  free(phi);
  free(SSHF);
  free(SS);
  free(GG);
  free(GG_new);
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_stls_arrays(input *in, double *xx, 
			    double *phi, double *SSHF, bool verbose){

  // Print on screen the parameter used to solve the STLS equation
  printf("------ Parameters used in the solution -------------\n");
  printf("Quantum degeneracy parameter: %f\n", in->Theta);
  printf("Quantum coupling parameter: %f\n", in->rs);
  printf("Chemical potential (low and high bound): %f %f\n", 
	 in->mu_lo, in->mu_hi);
  printf("Wave-vector cutoff: %f\n", in->xmax);
  printf("Wave-vector resolutions: %f\n", in->dx);
  printf("Number of Matsubara frequencies: %d\n", in->nl);
  printf("Maximum number of iterations: %d\n", in->nIter);
  printf("Error for convergence: %.5e\n", in->err_min_iter);
  printf("Number of threads: %d\n", omp_get_max_threads());
  printf("----------------------------------------------------\n");
 
  // Chemical potential
  if (in->Theta > 0) {
    if (verbose) printf("Chemical potential calculation: ");
    in->mu = compute_chemical_potential(*in);
    if (verbose) printf("Done. Chemical potential: %.8f\n", in->mu);
  }
  
  // Wave-vector grid
  if (verbose) printf("Wave-vector grid initialization: ");
  wave_vector_grid(xx, in);
  if (verbose) printf("Done.\n");
  
  // Normalized ideal Lindhard density response
  if (in->Theta > 0) {
    if (verbose) printf("Normalized ideal Lindhard density calculation:\n");
    compute_idr(phi, xx, *in, verbose);
    if (verbose) printf("Done.\n");
  }
  
  // Static structure factor in the Hartree-Fock approximation
  if (verbose) printf("Static structure factor in the Hartree-Fock approximation: ");
  if (in->Theta == 0) {
    compute_ssf_HF_zero_temperature(SSHF, xx, *in);
  }
  else {
    compute_ssf_HF_finite_temperature(SSHF, xx, *in);
  }
  if (verbose) printf("Done.\n");

}

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void wave_vector_grid(double *xx, input *in){
 
  xx[0] = 0.0;
  for (int ii=1; ii < in->nx; ii++) xx[ii] = xx[ii-1] + in->dx;

  // Ensure consistency between input and grid generated by the code
  in->xmax = xx[in->nx - 1];

}

// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A TWO-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx2(int xx, int yy, int x_size) {
  return (yy * x_size) + xx;
}


// -------------------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY AT FINITE TEMPERATURE
// -------------------------------------------------------------------------------------

struct idr_params {

  double xx;
  double mu;
  double Theta;
  double ll;

};

// Ideal density response for all matsubara frequencies specified in input
void compute_idr(double *phi, double *xx,  input in, bool verbose) {

  // Temporary array to store results
  double *phil = malloc( sizeof(double) * in.nx);
  
  // Loop over the Matsubara frequency
  for (int ll=0; ll<in.nl; ll++){
    // Compute lindhard density
    compute_idr_one_frequency(phil, xx, ll, in);
    // Fill output array
    for (int ii=0; ii<in.nx; ii++){
      phi[idx2(ii,ll,in.nx)] = phil[ii];
    }    
  }
  
  // Free memory
  free(phil);

}

// Ideal density response for one mastubara frequency
void compute_idr_one_frequency(double *phil, double *xx,  int ll, input in) {

  double err;
  size_t nevals;

  // Integration workspace 
  gsl_integration_cquad_workspace *wsp 
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  if (ll == 0) ff_int.function = &idr_partial_x0;
  else ff_int.function = &idr_partial_xl;

  // Normalized ideal Lindhard density 
  for (int ii = 0; ii < in.nx; ii++) {
    
    struct idr_params phixlp = {xx[ii], in.mu, in.Theta, ll};
    ff_int.params = &phixlp;
    gsl_integration_cquad(&ff_int, 
			  xx[0], xx[in.nx-1], 
			  0.0, QUAD_REL_ERR, 
			  wsp, 
			  &phil[ii], &err, &nevals);

  }

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

// Partial ideal density response (frequency = l, vector = x)
double idr_partial_xl(double yy, void *pp) {

  struct idr_params *params = (struct idr_params*)pp;
  double xx = (params->xx);
  double mu = (params->mu);
  double Theta = (params->Theta);
  double ll = (params->ll);
  double yy2 = yy*yy, xx2 = xx*xx, txy = 2*xx*yy, 
    tplT = 2*M_PI*ll*Theta, tplT2 = tplT*tplT;
  
  if (xx > 0.0) {
    return 1.0/(2*xx)*yy/(exp(yy2/Theta - mu) + 1.0)
      *log(((xx2+txy)*(xx2+txy) + tplT2)/((xx2-txy)*(xx2-txy) + tplT2));
  }
  else {
    return 0;
  }

}


// Partial ideal density response (frequency = 0, vector = x)
double idr_partial_x0(double yy, void *pp) {

  struct idr_params *params = (struct idr_params*)pp;
  double xx = (params->xx);
  double mu = (params->mu);
  double Theta = (params->Theta);
  double yy2 = yy*yy, xx2 = xx*xx, xy = xx*yy;

  if (xx > 0.0){

    if (xx < 2*yy){
      return 1.0/(Theta*xx)*((yy2 - xx2/4.0)*log((2*yy + xx)/(2*yy - xx)) + xy)
        *yy/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0);
    }
    else if (xx > 2*yy){
      return 1.0/(Theta*xx)*((yy2 - xx2/4.0)*log((2*yy + xx)/(xx - 2*yy)) + xy)
        *yy/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0);
    }
    else {
      return 1.0/(Theta)*yy2/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0);;
    }
  }

  else{
    return (2.0/Theta)*yy2/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0);
  }

}

// -------------------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY AT ZERO TEMPERATURE
// -------------------------------------------------------------------------------------


// Real part of the ideal density response
double idr_re_zero_temperature(double xx, double Omega) {

  double x_2 = xx/2.0;
  double Omega_2x = Omega/(2.0*xx);
  double sum_factor = x_2 + Omega_2x;
  double diff_factor = x_2 - Omega_2x;
  double sum_factor2 = sum_factor*sum_factor;
  double diff_factor2 = diff_factor*diff_factor;
  double log_sum_arg;
  double log_diff_arg;
  double adder1 = 0.0;
  double adder2 = 0.0;

  if (sum_factor != 1.0) {
    log_sum_arg = (sum_factor + 1.0)/(sum_factor - 1.0);
    if (log_sum_arg < 0.0) log_sum_arg = -log_sum_arg;
    adder1 = 1.0/(4.0*xx)*(1.0 - sum_factor2)*log(log_sum_arg);
  }

  if (diff_factor != 1.0 && diff_factor != -1.0) {
    log_diff_arg = (diff_factor + 1.0)/(diff_factor - 1.0);
    if (log_diff_arg < 0.0) log_diff_arg = -log_diff_arg;
    adder2 = 1.0/(4.0*xx)*(1.0 - diff_factor2)*log(log_diff_arg);
  }
  
  return 0.5 + adder1 + adder2;
  
}

// Imaginary part of the ideal density response
double idr_im_zero_temperature(double xx, double Omega) {

  double x_2 = xx/2.0;
  double Omega_2x = Omega/(2.0*xx);
  double sum_factor = x_2 + Omega_2x;
  double diff_factor = x_2 - Omega_2x;
  double sum_factor2 = sum_factor*sum_factor;
  double diff_factor2 = diff_factor*diff_factor;
  double adder1 = 0.0;
  double adder2 = 0.0;

  if (sum_factor2 < 1.0) {
    adder1 = 1 - sum_factor2;
  }

  if (diff_factor2 < 1.0) {
    adder2 = 1 - diff_factor2;
  }
  
  return -M_PI/(4.0*xx) * (adder1 - adder2);
  
}


// Frequency derivative of the real part of the ideal density response
double idrp_re_zero_temperature(double xx, double Omega) {

  double x_2 = xx/2.0;
  double Omega_2x = Omega/(2.0*xx);
  double sum_factor = x_2 + Omega_2x;
  double diff_factor = x_2 - Omega_2x;
  double log_sum_arg;
  double log_diff_arg;
  double adder1 = 0.0;
  double adder2 = 0.0;

  if (sum_factor != 1.0) {
    log_sum_arg = (sum_factor + 1.0)/(sum_factor - 1.0);
    if (log_sum_arg < 0.0) log_sum_arg = -log_sum_arg;
    adder1 = 1.0/(4.0*xx*xx)*(1.0 - sum_factor*log(log_sum_arg));
  }

  if (diff_factor != 1.0 && diff_factor != -1.0) {
    log_diff_arg = (diff_factor + 1.0)/(diff_factor - 1.0);
    if (log_diff_arg < 0.0) log_diff_arg = -log_diff_arg;
    adder2 = -1.0/(4.0*xx*xx)*(1.0 - diff_factor*log(log_diff_arg));
  }
  
  return adder1 + adder2;
  
}


// --------------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// --------------------------------------------------------------------------

// Static structure factor from the fluctuation-dissipation theorem
void compute_ssf_stls(double *SS, double *SSHF, double *GG, 
		      double *phi, double *xx, input in){

  if (in.Theta == 0) {
    compute_ssf_stls_zero_temperature(SS, SSHF, GG, xx, in);
  }
  else {
    compute_ssf_stls_finite_temperature(SS, SSHF, GG, phi, xx, in);
  }

}


// --------------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR AT FINITE TEMPERATURE
// --------------------------------------------------------------------------

// Static structure factor from the sum over the Matsubara frequencies
void compute_ssf_stls_finite_temperature(double *SS, double *SSHF, double *GG, 
					 double *phi, double *xx, input in){

  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double ff = 4*lambda*in.rs/M_PI;
  double xx2, BB, BB_tmp, BB_den, phixl;

  for (int ii=0; ii<in.nx; ii++){

    if (xx[ii] > 0.0){

      xx2 = xx[ii]*xx[ii];
      BB = 0.0;
      
      for (int ll=0; ll<in.nl; ll++){
	phixl = phi[idx2(ii,ll,in.nx)];
	BB_den = 1.0 + ff/xx2*(1 - GG[ii])*phixl;
	BB_tmp = phixl*phixl/BB_den;
	if (ll>0) BB_tmp *= 2.0;
	BB += BB_tmp;
	
      }
      
      SS[ii] = SSHF[ii]
	- 3.0/2.0*ff/xx2*in.Theta*(1- GG[ii])*BB;

    }
    else
      SS[ii] = 0.0;

  }

}

// Static structure factor within the Hartree-Fock approximation
struct ssf_HF_finite_temperature_params {

  double xx;
  double mu;
  double Theta;

};

void compute_ssf_HF_finite_temperature(double *SS,  double *xx,  input in){

  double err;
  size_t nevals;

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &ssf_HF_finite_temperature;

  // Static structure factor in the Hartree-Fock approximation
  for (int ii = 0; ii < in.nx; ii++) {


    struct ssf_HF_finite_temperature_params ssf_HF_p = {xx[ii], in.mu, in.Theta};
    ff_int.params = &ssf_HF_p;
    gsl_integration_cquad(&ff_int,
			  xx[0], xx[in.nx-1],
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &SS[ii], &err, &nevals);

    SS[ii] += 1.0;

  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

double ssf_HF_finite_temperature(double yy, void* pp) {

  struct idr_params *params = (struct idr_params*)pp;
  double xx = (params->xx);
  double mu = (params->mu);
  double Theta = (params->Theta);
  double yy2 = yy*yy, ypx = yy + xx, ymx = yy - xx;
 
  if (xx > 0.0){
    return -3.0*Theta/(4.0*xx)*yy/(exp(yy2/Theta - mu) + 1.0)
      *log((1 + exp(mu - ymx*ymx/Theta))/(1 + exp(mu - ypx*ypx/Theta)));
  }
  else {
    return -3.0*yy2/((1.0 + exp(yy2/Theta - mu))*(1.0 + exp(yy2/Theta - mu)));
  }

}



// ------------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR AT ZERO TEMPERATURE
// ------------------------------------------------------------------------

// Static structure factor from the integral over the frequencies
struct ssf_stls_zero_temperature_params {

  double xx;
  double rs;
  double GG;
  
};


void compute_ssf_stls_zero_temperature(double *SS, double *SSHF, double *GG, 
				       double *xx, input in){

  double ssf_wp;
  double err;
  double int_lo;
  double int_hi;
  size_t nevals;

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
  
  // Integration function
  gsl_function ff_int;
  ff_int.function = &ssf_stls_zero_temperature;
  
  for (int ii=0; ii<in.nx; ii++){

    if (xx[ii] > 0.0){

      // Integration limits
      if (xx[ii] < 2.0) int_lo = 0.0;
      else int_lo = xx[ii]*(xx[ii] - 2.0);
      int_hi = xx[ii]*(2.0 + xx[ii]);

      // Compute integral
      struct ssf_stls_zero_temperature_params ssf_p = {xx[ii], in.rs, GG[ii]};
      ff_int.params = &ssf_p;
      gsl_integration_cquad(&ff_int,
      			    int_lo, int_hi,
      			    0.0, QUAD_REL_ERR,
      			    wsp,
      			    &SS[ii], &err, &nevals);

      // Add plasmon contribution
      ssf_wp = ssf_plasmon(xx[ii], in);
      if (ssf_wp >= 0.0) SS[ii] += ssf_wp;

    }
    else {
      SS[ii] = 0.0;
    }

    // Add Hartree-Fock contribution
    SS[ii] += SSHF[ii];
    
  }

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}


double ssf_stls_zero_temperature(double Omega, void* pp) {

  struct ssf_stls_zero_temperature_params *params = (struct ssf_stls_zero_temperature_params*)pp;
  double xx = (params->xx);
  double rs = (params->rs);
  double GG = (params->GG);
  double xx2 = xx*xx;
  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0); 
  double ff = (4.0*lambda*rs)/(M_PI*xx2);
  double phi0_im;
  double phi0_re;
  double fact_re, fact_re2;
  double fact_im, fact_im2;

  phi0_re = idr_re_zero_temperature(xx, Omega);
  phi0_im = idr_im_zero_temperature(xx, Omega);

  fact_re = 1 + ff*(1 - GG)*phi0_re;
  fact_im = ff*(1 - GG)*phi0_im;
  fact_re2 = fact_re*fact_re;
  fact_im2 = fact_im*fact_im;

  return 3.0/(2.0*M_PI)*phi0_im*(1.0/(fact_re2 + fact_im2) - 1.0);

}

// Static structure factor within the Hartree-Fock approximation
void compute_ssf_HF_zero_temperature(double *SS,  double *xx,  input in){

  for (int ii = 0; ii < in.nx; ii++) {
    if (xx[ii] < 2.0) {
      SS[ii] = (1.0/16.0)*xx[ii]*(12.0 - xx[ii]*xx[ii]);
    }
    else {
      SS[ii] = 1.0;
    }
  }
  
}


// Plasmon contribution to the static structure factor
struct eps_params {

  double xx;
  double rs;

};


double ssf_plasmon(double xx, input in) {
  
  // Variables
  double w_co = xx*xx + 2*xx;
  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double ff = (4.0*lambda*in.rs)/(M_PI*xx*xx);
  double dx;
  double w_lo, w_min, w_hi;
  double eps_lo, eps_min, eps_hi;
  int status, iter;

  // Set-up function
  gsl_function ff_minimizer;
  ff_minimizer.function = &dr_mod_zero_temperature;
  struct eps_params epsp = {xx, in.rs};
  ff_minimizer.params = &epsp;
    
  // Get approximate location of the mimimum (to initialize the minimizer)
  dx = w_co;
  w_lo = -dx;
  w_min = 0.0;
  w_hi = dx;
  eps_lo = 0.0;
  eps_min = 1.0;
  eps_hi = 0.0;
  while (eps_min > eps_lo || eps_min > eps_hi) {
    w_lo = w_min;
    w_min = w_hi;
    w_hi += dx;
    eps_lo = dr_mod_zero_temperature(w_lo, &epsp);
    eps_min = dr_mod_zero_temperature(w_min, &epsp);
    eps_hi = dr_mod_zero_temperature(w_hi, &epsp); 
  }
  
  // Set-up minimizer
  const gsl_min_fminimizer_type * minit = gsl_min_fminimizer_goldensection;
  gsl_min_fminimizer * mini = gsl_min_fminimizer_alloc(minit);
  gsl_min_fminimizer_set(mini, &ff_minimizer, w_min, w_lo, w_hi);

  // Solve dispersion relation to find the plasmon frequency
  iter = 0;
  do
  {
    
    // Solver iteration
    status = gsl_min_fminimizer_iterate(mini);

    // Get solver status
    w_min = gsl_min_fminimizer_minimum(mini);
    w_lo = gsl_min_fminimizer_x_lower(mini);
    w_hi = gsl_min_fminimizer_x_upper(mini);
    status = gsl_min_test_interval(w_lo, w_hi, ROOTMIN_ABS_ERR, 0.0);

    // Update iteration counter
    iter++;

  }
  while (status == GSL_CONTINUE && iter < ROOTMIN_MAX_ITER);

  // Free memory
  gsl_min_fminimizer_free(mini);

  // Dielectric response at the minimum
  eps_min = dr_mod_zero_temperature(w_min, &epsp);
  
  // Check value of the function at the minima and return accordingly
  /* if (status == GSL_SUCCESS && eps_min < ROOTMIN_REL_ERR) { */
  /*   printf("--------------------------------\n"); */
  /*   printf("%f %f %.10f %.10f\n", xx, w_min, eps_min, 1.5/(ff*fabs(drp_re_zero_temperature(xx, w_min, in.rs)))); */
  /* } */
  if (status == GSL_SUCCESS && eps_min < ROOTMIN_ABS_ERR) {
    return 1.5 / (ff*fabs(drp_re_zero_temperature(xx, w_min, in.rs)));
  } else {
    return -1; // No valid root was found
  }
  return 0;
  
}

// Real part of the dielectric response
double dr_re_zero_temperature(double xx, double Omega, double rs){

  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0); 
  double ff = (4.0*lambda*rs)/(M_PI*xx*xx);

  return 1 + ff*idr_re_zero_temperature(xx,  Omega);
    
}

// Imaginary part of the dielectric response
double dr_im_zero_temperature(double xx, double Omega, double rs){

  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0); 
  double ff = (4.0*lambda*rs)/(M_PI*xx*xx);

  return ff*idr_im_zero_temperature(xx,  Omega);
    
}

// Modulus of the dielectric response
double dr_mod_zero_temperature(double Omega, void *pp){

  struct eps_params *params = (struct eps_params*)pp;
  double xx = params->xx;
  double rs = params->rs;

  double eps_re = dr_re_zero_temperature(xx, Omega, rs); 
  double eps_im = dr_im_zero_temperature(xx, Omega, rs);
  
  return sqrt(eps_re*eps_re + eps_im*eps_im);
  
}

// Frequency derivative of the real part of the dielectric response
double drp_re_zero_temperature(double xx, double Omega, double rs){

  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0); 
  double ff = (4.0*lambda*rs)/(M_PI*xx*xx);

  return ff*idrp_re_zero_temperature(xx, Omega);
    
}
 
// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

struct slfc_params {

  double xx;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};


void compute_slfc(double *GG, double *SS, double *xx, input in) {

  double err;
  size_t nevals;
  
  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);    

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &slfc;

  // Static local field correction
  for (int ii = 0; ii < in.nx; ii++) {
    struct slfc_params slfcp = {xx[ii], ssf_sp_ptr, ssf_acc_ptr};
    ff_int.params = &slfcp;
    gsl_integration_cquad(&ff_int,
			  xx[0], xx[in.nx-1],
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &GG[ii], &err, &nevals);
  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(ssf_sp_ptr);
  gsl_interp_accel_free(ssf_acc_ptr);
  
}

double slfc(double yy, void* pp) {

  struct slfc_params* params = (struct slfc_params*)pp;
  double xx = (params->xx);
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  double yy2 = yy * yy, xx2 = xx * xx;

  if (xx > 0.0 && yy > 0.0){

    if (xx > yy){
      return -(3.0/4.0)* yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0)
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(xx - yy)));
    }
    else if (xx < yy) {
      return -(3.0/4.0) * yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0)
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(yy - xx)));
    }
    else {
      return -(3.0/4.0) * yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0);
    }

  }
  else
    return 0;



}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

struct uex_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};

double compute_internal_energy(double *SS, double *xx,  input in) {

  double err;
  size_t neval;
  double ie;
  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);  
  
  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_linear, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &uex;

  // Internal energy
  struct uex_params uexp = {ssf_sp_ptr, ssf_acc_ptr};
  ff_int.params = &uexp;  
  gsl_integration_cquad(&ff_int,
			xx[0], xx[in.nx-1],
			0.0, QUAD_REL_ERR,
			wsp,
			&ie, &err, &neval);
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(ssf_sp_ptr);
  gsl_interp_accel_free(ssf_acc_ptr);

  // Output
  return ie/(M_PI*in.rs*lambda);

}

double uex(double yy, void* pp) {

  struct uex_params* params = (struct uex_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);

  return gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1;
}


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------


// write text files for output
void write_text_stls(double *SS, double *GG, double *phi, 
		       double *SSHF, double *xx, input in){


    FILE* fid;
    
    // Output for SSF
    char out_name[100];
    sprintf(out_name, "ssf_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static structure factor\n");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nx; ii++)
        fprintf(fid, "%.8e %.8e\n", xx[ii], SS[ii]);

    fclose(fid);

    // Output for SLFC
    sprintf(out_name, "slfc_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static local field correction\n");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nx; ii++)
        fprintf(fid, "%.8e %.8e\n", xx[ii], GG[ii]);

    fclose(fid);

    // Output for static density response
    sprintf(out_name, "sdr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
      perror("Error while creating the output file for the static density response\n");
      exit(EXIT_FAILURE);
    }
    double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
    double ff = 4*lambda*in.rs/M_PI;
    double sdr;
    for (int ii=0 ; ii<in.nx; ii++){
	sdr = -(3.0/2.0)*in.Theta*phi[idx2(ii,0,in.nx)]/
	  (1.0 + ff/(xx[ii]*xx[ii])*(1.0 - GG[ii])*phi[idx2(ii,0,in.nx)]);
	fprintf(fid, "%.8e %.8e\n", xx[ii], sdr);
      }
    fclose(fid);

    // Output for ideal Lindhard density response
    sprintf(out_name, "idr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
      perror("Error while creating the output file for the ideal density response");
      exit(EXIT_FAILURE);
    }
    for (int ii=0; ii<in.nx; ii++){
      for (int jj=0; jj<in.nl; jj++){
        fprintf(fid, "%.8e ", phi[idx2(ii,jj,in.nx)]);
      }
      fprintf(fid,"\n");
    }
    fclose(fid);

    // Output for static structure factor in the Hartree-Fock approximation
    sprintf(out_name, "ssfHF_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static structure factor (HF)");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nx; ii++)
        fprintf(fid, "%.8e %.8e\n", xx[ii], SSHF[ii]);

    fclose(fid);

    // Output for the interaction energy
    sprintf(out_name, "uint_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the interaction energy\n");
        exit(EXIT_FAILURE);
    }
    fprintf(fid, "%.8e\n", compute_internal_energy(SS, xx, in));
    fclose(fid);

}


// write binary file to use as initial guess (or restart)
void write_guess_stls(double *SS, double *GG, input in){

  // Name of output file
  char out_name[100];
  sprintf(out_name, "restart_rs%.3f_theta%.3f_%s.bin", in.rs, in.Theta, in.theory);

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(out_name, "wb");
  if (fid == NULL) {
    fprintf(stderr,"Error while creating file for initial guess or restart\n");
    exit(EXIT_FAILURE);
  }

  // Input data
  fwrite(&in.nx, sizeof(int), 1, fid);
  fwrite(&in.dx, sizeof(double), 1, fid);
  fwrite(&in.xmax, sizeof(double), 1, fid);

  // Static structure factor 
  fwrite(SS, sizeof(double), in.nx, fid);

  // Static local field correction
  fwrite(GG, sizeof(double), in.nx, fid);

  // Close binary file
  fclose(fid);

}


// read binary file to use as initial guess (or restart)
void read_guess_stls(double *SS, double *GG, input in){

  // Variables
  size_t it_read;
  int nx_file;
  double dx_file;
  double xmax_file;

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(in.stls_guess_file, "rb");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file for initial guess or restart\n");
    exit(EXIT_FAILURE);
  }

  // Initialize number of items read from input file
  it_read = 0;

  // Check that the data for the guess file is consistent
  it_read += fread(&nx_file, sizeof(int), 1, fid);
  it_read += fread(&dx_file, sizeof(double), 1, fid);
  it_read += fread(&xmax_file, sizeof(double), 1, fid);
  check_guess_stls(nx_file, dx_file, xmax_file, in, it_read, 3,
		   fid, true, true, false);   
  
  // Static structure factor in the Hartree-Fock approximation
  it_read += fread(SS, sizeof(double), nx_file, fid);
  
  // Static structure factor in the Hartree-Fock approximation
  it_read += fread(GG, sizeof(double), nx_file, fid);

  // Check that all items where read and the end-of-file was reached
  check_guess_stls(nx_file, dx_file, xmax_file, in, it_read,
		   2*nx_file + 3, fid, false, true, true);   
 
  // Close binary file
  fclose(fid);
	    
}


// Check consistency of the guess data
void check_guess_stls(int nx, double dx, double xmax, input in,
		      size_t it_read, size_t it_expected, FILE *fid,
		      bool check_grid, bool check_items, bool check_eof){

  int buffer;
  double tol = 1e-10;
  
  // Check that the grid in the guess data is consistent with input
  if (check_grid) {
    
    if (nx != in.nx || (dx-in.dx) > tol || (xmax-in.xmax) > tol){
      fprintf(stderr,"Grid from guess file is incompatible with input\n");
      fprintf(stderr,"Grid points (nx) : %d (input), %d (file)\n", in.nx, nx);
      fprintf(stderr,"Resolution (dx)  : %.16f (input), %.16f (file)\n", in.dx, dx);
      fprintf(stderr,"Cutoff (xmax)    : %.16f (input), %.16f (file)\n", in.xmax, xmax);
      fclose(fid);
      exit(EXIT_FAILURE);
    }

  }

  // Check that all the expected items where read
  if (check_items) {
    if (it_read != it_expected ) {
      fprintf(stderr,"Error while reading file for initial guess or restart.\n");
      fprintf(stderr,"%ld Elements expected, %ld elements read\n", it_read, it_expected);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
  }
  
  // Check for end of file
  if (check_eof){
    it_read = fread(&buffer, sizeof(int), 1, fid); // Trigger end-of-file activation
    if (!feof(fid)) {
      fprintf(stderr,"Error while reading file for initial guess or restart.\n");
      fprintf(stderr,"Expected end of file, but there is still data left to read.\n");
      fclose(fid);
      exit(EXIT_FAILURE);
    }
  }
  
  
}

