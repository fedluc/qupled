#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "solvers.h"
#include "utils.h"
#include "chemical_potential.h"
#include "stls.h"
#include "stls_iet.h"
#include "qstls.h"
#include "dynamic_stls.h"
#include "dynamic_qstls.h"
#include "dynamic_qstls_iet.h"

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC PROPERTIES OF QSTLS SCHEME
// -------------------------------------------------------------------

void compute_dynamic_qstls_iet(input in, bool verbose) {

  // Arrays
  double *WW = NULL;
  double *phi_re = NULL;
  double *phi_im = NULL;
  double *psi_re = NULL;
  double *psi_im = NULL;
  double *SSn = NULL;
  double *SS = NULL;
  double *xx = NULL;
  double *bf = NULL;
  
  // Safeguard
  if (in.Theta == 0) {
    printf("Ground state calculations of the dynamic properties"
	   " are not yet implemented.");
    exit(EXIT_FAILURE);
  }
      
  // Get the size of the frequency grid
  get_frequency_grid_size(&in);

  // Static structure factor (this also defines the wave-vector grid)
  if (verbose) printf("Static structure factor (from file): ");
  get_ssf(&SS, &xx, &in);
  if (verbose) printf("Done.\n");

  // Allocate arrays
  alloc_dynamic_qstls_iet_arrays(in, &WW, &phi_re, &phi_im, &psi_re,
				 &psi_im, &SSn);
  
  // Chemical potential and frequency grid
  init_fixed_dynamic_stls_arrays(&in, WW, verbose);
  
  // Ideal density response
  if (verbose) printf("Normalized ideal Lindhard density calculation: ");
  compute_dynamic_idr_iet(phi_re, phi_im, WW, xx, in);
  if (verbose) printf("Done.\n");
  
  // Bridge function term
  if (verbose) printf("Bridge function Fourier transform: ");
  get_bf(&bf, xx, in);
  if (verbose) printf("Done.\n");
  
  // Auxiliary density response
  if (verbose) printf("Auxiliary density calculation: ");
  compute_dynamic_adr_iet(psi_re, psi_im, phi_re, phi_im, WW, SS,
			  bf, xx, in);
  if (verbose) printf("Done.\n");

  /* // Dynamic structure factor */
  /* if (verbose) printf("Dynamic structure factor calculation: "); */
  /* compute_dsf_qstls(SSn, phi_re, phi_im, psi_re, psi_im, WW, in); */
  /* if (verbose) printf("Done.\n"); */
  
  /* // Output to file */
  /* if (verbose) printf("Writing output files: "); */
  /* write_text_dynamic_qstls(SSn, WW, psi_re, psi_im, in); */
  /* if (verbose) printf("Done.\n"); */

  // Free memory
  free_dynamic_stls_arrays(WW, phi_re, phi_im, SSn);
  free_dynamic_qstls_arrays(psi_re, psi_im, SS, xx);
  free(bf);
  
 
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_qstls_iet_arrays(input in, double **WW,
				    double **phi_re, double **phi_im,
				    double **psi_re, double **psi_im,
				    double **SSn){

  *WW = malloc( sizeof(double) * in.nW);
  if (*WW == NULL) {
    fprintf(stderr, "Failed to allocate memory for the frequency grid\n");
    exit(EXIT_FAILURE);
  }

  *SSn = malloc( sizeof(double) * in.nW);
  if (*SSn == NULL) {
    fprintf(stderr, "Failed to allocate memory for the dynamic structure factor\n");
    exit(EXIT_FAILURE);
  }
  
  *phi_re = malloc( sizeof(double) * in.nx * in.nW);
  if (*phi_re == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *phi_im = malloc( sizeof(double) * in.nx * in.nW);
  if (*phi_im == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *psi_re = malloc( sizeof(double) * in.nx * in.nW);
  if (*psi_re == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the auxiliary density response\n");
    exit(EXIT_FAILURE);
  }
  
  *psi_im = malloc( sizeof(double) * in.nx * in.nW);
  if (*psi_im == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the auxiliary density response\n");
    exit(EXIT_FAILURE);
  }
  
}


// -------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE BRIDGE FUNCTION
// -------------------------------------------------------------------

void get_bf(double **bf, double *xx, input in){

  // Allocate array to store the bridge function
  *bf = malloc( sizeof(double) * in.nx);
  if (*bf == NULL) {
    fprintf(stderr, "Failed to allocate memory for the bridge function\n");
    exit(EXIT_FAILURE);
  }

  // Compute bridge function
  compute_bridge_function(*bf, xx, in);
  
}


// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE IDEAL DENSITY RESPONSE
// ------------------------------------------------------------------

void compute_dynamic_idr_iet(double *phi_re, double *phi_im,
			     double *WW, double *xx, input in) {

  // Allocate temporary arrays
  double *phi_re_tmp = malloc( sizeof(double) * in.nW);
  double *phi_im_tmp = malloc( sizeof(double) * in.nW);
  if (phi_re_tmp == NULL ||
      phi_im_tmp == NULL){
    fprintf(stderr, "Failed to allocate memory for calculation"
	    " of the ideal density response function\n");
    exit(EXIT_FAILURE);
  }

  // Loop over the wave-vector grid
  for (int ii=0; ii<in.nx; ii++) {

    // Ideal density response for one wave-vector
    in.dyn_xtarget = xx[ii];
    compute_dynamic_idr(phi_re_tmp, phi_im_tmp, WW, in);

    // Copy temporary arrays
    for (int jj=0; jj<in.nW; jj++) {
      phi_re[idx2(ii,jj,in.nx)] = phi_re_tmp[jj];
      phi_im[idx2(ii,jj,in.nx)] = phi_im_tmp[jj];
    }
    
  }
  
  // Free memory
  free(phi_re_tmp);
  free(phi_im_tmp);

}


// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE AUXILIARY DENSITY RESPONSE
// ------------------------------------------------------------------

// Ideal density response (real and imaginary part)
void compute_dynamic_adr_iet(double *psi_re, double *psi_im,
			     double *phi_re, double *phi_im,
			     double *WW, double *SS,
			     double *bf, double *xx,
			     input in) {

  // Real component
  compute_dynamic_adr_iet_re(psi_re, phi_re, WW, SS, bf, xx, in);
  
  // Imaginary component
  //compute_dynamic_adr_im_lev1(psi_im, WW, SS, xx, in);
  
}


// ------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE REAL PART OF THE AUXILIARY
// DENSITY RESPONSE
// ------------------------------------------------------------------

// Real part of the auxiliary density response (iterations)
void compute_dynamic_adr_iet_re(double *psi_re, double *phi_re,
				double *WW, double *SS,
				double *bf, double *xx,
				input in){

  int iter_counter;
  double iter_err;
  
  // Allocate temporary arrays
  double *psi_re_new = malloc( sizeof(double) * in.nx * in.nW);
  double *psi_re_fixed = malloc( sizeof(double) * in.nx * in.nW);
  if (psi_re == NULL){
    fprintf(stderr, "Failed to allocate memory for calculation"
	    " of the real part of the auxiliary density"
	    " response function\n");
    exit(EXIT_FAILURE);
  }

  // Initialize the auxiliary density response
  for (int ii=0; ii<in.nx; ii++){
    for (int jj=0; jj<in.nW; jj++){
      psi_re[idx2(ii,jj,in.nx)] = 0.0;
      psi_re_fixed[idx2(ii,jj,in.nx)] = INFINITY;
    }
  }

  // Iterations
  iter_counter = 0;
  iter_err = 1;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter) {

    // Update auxiliary density response
    compute_dynamic_adr_iet_re_lev1(psi_re_new, psi_re, psi_re_fixed,
				    phi_re, WW, SS, bf, xx, in);

    // Update diagnostics
    iter_counter++;
    iter_err = adr_iet_err(psi_re, psi_re_new, in);
    adr_iet_update(psi_re, psi_re_new, in);

    if (true) {
      printf("--- iteration %d ---\n", iter_counter);
      printf("Residual error: %.5e\n", iter_err);
      fflush(stdout);
    }
    
  }
  
  // Free memory
  free(psi_re_fixed);
  free(psi_re_new);

}


// Relative error for the iterations
double adr_iet_err(double *psi_re, double *psi_re_new, input in){

  double err = 0.0;
  double err_tmp;
  
  for (int ii=0; ii<in.nx; ii++){
    err_tmp = psi_re[idx2(ii,0,in.nx)] - psi_re_new[idx2(ii,0,in.nx)];
    err += err_tmp*err_tmp;
  }

  return sqrt(err);
  
}

// Update iterative solution
void adr_iet_update(double *psi_re, double *psi_re_new, input in){

  for (int ii=0; ii<in.nx; ii++){
    psi_re[idx2(ii,0,in.nx)] = in.a_mix * psi_re_new[idx2(ii,0,in.nx)]
      + (1 - in.a_mix)*psi_re[idx2(ii,0,in.nx)];
  }

}


// Real part of the auxiliary density response (level 1)

struct adr_re_lev1_params {

  gsl_spline *int_lev1_1_sp_ptr;
  gsl_interp_accel *int_lev1_1_acc_ptr;
  gsl_spline *int_lev1_2_sp_ptr;
  gsl_interp_accel *int_lev1_2_acc_ptr;

};

void compute_dynamic_adr_iet_re_lev1(double *psi_re_new, double *psi_re,
				     double *psi_re_fixed, double *phi_re,
				     double *WW, double *SS,
				     double *bf, double *xx,
				     input in) {

  // Parallel calculations
  #pragma omp parallel
  {

    bool compute_fixed = true;
    double err;
    size_t nevals;
    double *int_lev1_1  = malloc( sizeof(double) * in.nx);
    double *int_lev1_2  = malloc( sizeof(double) * in.nx);
    if (int_lev1_1 == NULL ||
	int_lev1_2 == NULL){
      fprintf(stderr, "Failed to allocate memory for calculation"
	      " of the real part of the auxiliary density"
	      " response function\n");
      exit(EXIT_FAILURE);
    }
    
    // Declare accelerator and spline objects
    gsl_spline *int_lev1_1_sp_ptr;
    gsl_interp_accel *int_lev1_1_acc_ptr;
    gsl_spline *int_lev1_2_sp_ptr;
    gsl_interp_accel *int_lev1_2_acc_ptr;
    
    // Allocate the accelerator and the spline objects
    int_lev1_1_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    int_lev1_1_acc_ptr = gsl_interp_accel_alloc();
    int_lev1_2_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    int_lev1_2_acc_ptr = gsl_interp_accel_alloc();
    
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);

    if (psi_re_fixed[0] != INFINITY) compute_fixed = false;
    
    // Loop over the wave-vectors
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){

      // Loop over the frequencies
      for (int jj=0; jj<in.nW; jj++) {

	// Integration function
	gsl_function ff_int_lev1;
	ff_int_lev1.function = &adr_iet_re_lev1_partial_xW;

	// Integrand (part 1)
	compute_dynamic_adr_iet_re_lev1_1(int_lev1_1, psi_re, phi_re, SS, bf, in);
	gsl_spline_init(int_lev1_1_sp_ptr, xx, int_lev1_1, in.nx);
	  
	// Integrand (part 2)
	if (compute_fixed) {
	  compute_dynamic_adr_iet_re_lev2(int_lev1_2, WW[jj], xx[ii], SS, xx, in);
	  compute_dynamic_adr_iet_re_lev1_2(int_lev1_2, psi_re_fixed, jj, in, false);
	}
	else {
	  compute_dynamic_adr_iet_re_lev1_2(int_lev1_2, psi_re_fixed, jj, in, true);
	}
	gsl_spline_init(int_lev1_2_sp_ptr, xx, int_lev1_2, in.nx);
	      
	// Integral over w
	struct adr_re_lev1_params plev1 = {int_lev1_1_sp_ptr,
					   int_lev1_1_acc_ptr,
					   int_lev1_2_sp_ptr,
					   int_lev1_2_acc_ptr};
	ff_int_lev1.params = &plev1;
	gsl_integration_cquad(&ff_int_lev1,
			      xx[0], xx[in.nx-1],
			      0.0, QUAD_REL_ERR,
			      wsp,
			      &psi_re_new[idx2(ii,jj,in.nx)],
			      &err, &nevals);
      }
    }
    
    // Free memory
    free(int_lev1_1);
    free(int_lev1_2);
    gsl_integration_cquad_workspace_free(wsp);
    gsl_spline_free(int_lev1_1_sp_ptr);
    gsl_interp_accel_free(int_lev1_1_acc_ptr);
    gsl_spline_free(int_lev1_2_sp_ptr);
    gsl_interp_accel_free(int_lev1_2_acc_ptr);
    
  }
  
}

// Integrand for level 1 of the real auxiliary density response (part 1)
void compute_dynamic_adr_iet_re_lev1_1(double *int_lev1_1, double *psi_re,
				       double *phi_re, double *SS,
				       double *bf, input in){

  double psi_phi;
  
  for (int ii=0; ii<in.nx; ii++){
    if (ii==0) int_lev1_1[ii] = 0.0;
    else {
      psi_phi = psi_re[idx2(ii,0,in.nx)]/phi_re[idx2(ii,0,in.nx)];
      int_lev1_1[ii] = SS[ii]*(1.0 - bf[ii]) - psi_phi*(SS[ii] - 1.0);
    }
  }
  
}

// Integrand for level 1 of the real auxiliary density response (part 2)
void compute_dynamic_adr_iet_re_lev1_2(double *int_lev1_2, double *psi_re_fixed,
				       int jj, input in, bool read){

  if (read)
    for (int ii=0; ii<in.nx; ii++)
      int_lev1_2[ii] = psi_re_fixed[idx2(ii,jj,in.nx)];
  else 
    for (int ii=0; ii<in.nx; ii++)
      psi_re_fixed[idx2(ii,jj,in.nx)] = int_lev1_2[ii];
  
}


// Integrand for level 1 of the real auxiliary density response (vector = x, frequency = W)
double adr_iet_re_lev1_partial_xW(double ww, void* pp) {
  
  struct adr_re_lev1_params* params = (struct adr_re_lev1_params*)pp;
  gsl_spline* int_lev1_1_sp_ptr = (params->int_lev1_1_sp_ptr);
  gsl_interp_accel* int_lev1_1_acc_ptr = (params->int_lev1_1_acc_ptr);
  gsl_spline* int_lev1_2_sp_ptr = (params->int_lev1_2_sp_ptr);
  gsl_interp_accel* int_lev1_2_acc_ptr = (params->int_lev1_2_acc_ptr);
  double fflv1_1 = gsl_spline_eval(int_lev1_1_sp_ptr, ww, int_lev1_1_acc_ptr);
  double fflv1_2 = gsl_spline_eval(int_lev1_2_sp_ptr, ww, int_lev1_2_acc_ptr);

  if (ww == 0.0)
    return 0;
  else
    return fflv1_1*fflv1_2/ww;

}


// Real part of the auxiliary density response (level 2)

struct adr_re_lev2_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *int_lev2_sp_ptr;
  gsl_interp_accel *int_lev2_acc_ptr;
  
};

void compute_dynamic_adr_iet_re_lev2(double *int_lev1, double WW,
				     double xx, double *SS,
				     double *ww, input in) {

  double err;
  size_t nevals;
  double u_min;
  double u_max;
  double w_max = ww[in.nx-2];
  double *int_lev2  = malloc( sizeof(double) * in.nx);
  if (int_lev2 == NULL){
    fprintf(stderr, "Failed to allocate memory for calculation"
	    " of the real part of the auxiliary density"
	    " response function\n");
    exit(EXIT_FAILURE);
  }
  
  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *int_lev2_sp_ptr;
  gsl_interp_accel *int_lev2_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  int_lev2_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  int_lev2_acc_ptr = gsl_interp_accel_alloc();
  
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
     
  // Integration function
  gsl_function ff_int_lev2;
  ff_int_lev2.function = &adr_iet_re_lev2_partial_xwW;

  // Interpolation of the static structure factor
  gsl_spline_init(ssf_sp_ptr, ww, SS, in.nx);
  
  // Loop over w (wave-vector)
  for (int ii=0; ii<in.nx; ii++) {

    // Integration limits
    u_min = ww[ii] - xx;
    if (u_min < 0.0) u_min = -u_min;
    u_max = ww[ii] + xx;
    if (u_max > w_max) u_max = w_max;
    
    // Inner integral
    compute_dynamic_adr_iet_re_lev3(int_lev2, WW, xx, ww[ii], ww, in);

    // Construct integrand
    gsl_spline_init(int_lev2_sp_ptr, ww, int_lev2, in.nx);
    
    // Integration over u (wave-vector squared)
    struct adr_re_lev2_params plev2 = {ssf_sp_ptr,
				       ssf_acc_ptr,
				       int_lev2_sp_ptr,
				       int_lev2_acc_ptr};
    
    ff_int_lev2.params = &plev2;
    gsl_integration_cquad(&ff_int_lev2,
			  u_min, u_max,
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &int_lev1[ii],
			  &err, &nevals);
  }
 
    

  // Free memory
  free(int_lev2);
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(ssf_sp_ptr);
  gsl_interp_accel_free(ssf_acc_ptr);
  gsl_spline_free(int_lev2_sp_ptr);
  gsl_interp_accel_free(int_lev2_acc_ptr);
  
}


// Integrand for level 2 of the real auxiliary density response (vectors = {x,w}, frequency = W)
double adr_iet_re_lev2_partial_xwW(double uu, void* pp) {
  
  struct adr_re_lev2_params* params = (struct adr_re_lev2_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  gsl_spline* int_lev2_sp_ptr = (params->int_lev2_sp_ptr);
  gsl_interp_accel* int_lev2_acc_ptr = (params->int_lev2_acc_ptr);
  double ssfm1 = gsl_spline_eval(ssf_sp_ptr, uu, ssf_acc_ptr) - 1.0;
  double fflv2 = gsl_spline_eval(int_lev2_sp_ptr, uu, int_lev2_acc_ptr);

  return uu*ssfm1*fflv2;
  
}

// Real part of the auxiliary density response (level 3)

struct adr_re_lev3_params {

  double mu;
  double Theta;
  double xx;
  double ww;
  double uu;
  double WW;
  
};

void compute_dynamic_adr_iet_re_lev3(double *int_lev2, double WW,
				     double xx, double ww,
				     double *uu, input in) {

  double err;
  size_t nevals;
 
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
    
  // Integration function
  gsl_function ff_int_lev3;
  if (WW == 0.0)
    ff_int_lev3.function = &adr_re_lev3_partial_xwu0;
  else
    ff_int_lev3.function = &adr_re_lev3_partial_xwuW;
  
  // Loop over u (wave-vector)
  for (int ii=0; ii<in.nx; ii++){

    // Integrate over q (wave-vector)
    struct adr_re_lev3_params plev3 = {in.mu,in.Theta, xx, ww,
				       uu[ii], WW};
    ff_int_lev3.params = &plev3;
    gsl_integration_cquad(&ff_int_lev3,
			  uu[0], uu[in.nx-1],
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &int_lev2[ii],
			  &err, &nevals);
  }

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

// Integrand for level 3 of the real  auxiliary density response (vectors = {x,w,u}, frequency = W)
double adr_iet_re_lev3_partial_xwuW(double qq, void* pp) {
  
  struct adr_re_lev3_params* params = (struct adr_re_lev3_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  double xx = (params->xx);
  double ww = (params->ww);
  double uu = (params->uu);
  double WW = (params->WW);
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double qq2 = qq*qq;
  double uu2 = uu*uu;
  double WW2 = WW*WW;
  double f1 = xx2 + ww2 - uu2 + 4.0*xx*qq;
  double f2 = xx2 + ww2 - uu2 - 4.0*xx*qq;
  double logarg = (f1*f1 - 4*WW2)/(f2*f2 - 4.0*WW2);

  if (logarg < 0) logarg = -logarg;
  
  return -(3.0/8.0)*qq/(exp(qq2/Theta - mu) + 1.0)*
    log(logarg);

}


// Integrand for level 3 of the real  auxiliary density response (vectors = {x,w,u}, frequency = 0)
double adr_iet_re_lev3_partial_xwu0(double qq, void* pp) {
  
  struct adr_re_lev3_params* params = (struct adr_re_lev3_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  double xx = (params->xx);
  double ww = (params->ww);
  double uu = (params->uu);
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double qq2 = qq*qq;
  double uu2 = uu*uu;
  double x2u2w2 = xx2 + ww2 - uu2;
  double f1 = x2u2w2 + 4.0*xx*qq;
  double f2 = x2u2w2 - 4.0*xx*qq;
  double logarg = (f1*f1)/(f2*f2);

  if (xx == 0 || qq == 0){
    return 0;
  }
  else {
    
    if (logarg < 0.0) logarg = -logarg;
    return  -(3.0/(4.0*Theta))
      *qq/(exp(qq2/Theta - mu)+ exp(-qq2/Theta + mu) + 2.0)
      *((qq2 - x2u2w2*x2u2w2/(16.0*xx2))*log(logarg)
	+ (qq/xx)*x2u2w2/2.0);
				     
  }
  

}

/* // --------------------------------------------------------------------- */
/* // FUNCTION USED TO COMPUTE THE DYNAMIC STRUCTURE FACTOR */
/* // --------------------------------------------------------------------- */

/* void compute_dsf_qstls(double *SSn, double *phi_re, double *phi_im, */
/* 		       double *psi_re, double *psi_im, */
/* 		       double *WW, input in){ */
    
/*   double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0); */
/*   double xx = in.dyn_xtarget; */
/*   double ff1 = 4.0*lambda*in.rs/(M_PI*xx*xx); */
/*   double ff2; */
/*   double numer; */
/*   double denom; */
/*   double denom_re; */
/*   double denom_im; */
  
/*   for (int ii=0; ii<in.nW; ii++){ */

/*     if (WW[ii] == 0.0) { */

/*       ff2 = in.Theta/(4.0*xx); */
/*       numer = (1.0 - ff1*psi_re[ii]) */
/* 	/(exp(xx*xx/(4.0*in.Theta) - in.mu) + 1) */
/* 	- 3.0/(4.0*xx)*ff1*phi_re[ii]*psi_im[ii]; */
/*       numer *= ff2; */
/*       denom_re = 1.0 + ff1 * (phi_re[ii] - psi_re[ii]); */
/*       denom = denom_re * denom_re; */
	
/*     } */
/*     else { */
      
/*       ff2 = 1.0/(1.0 - exp(-WW[ii]/in.Theta)); */
/*       numer = phi_im[ii] + ff1*(phi_re[ii]*psi_im[ii] - */
/* 				phi_im[ii]*psi_re[ii]); */
/*       numer *= (ff2/M_PI); */
/*       denom_re = 1.0 + ff1 * (phi_re[ii] - psi_re[ii]); */
/*       denom_im = ff1 * (phi_im[ii] - psi_im[ii]); */
/*       denom = denom_re*denom_re + denom_im*denom_im; */
      
/*     } */
    
/*     if (xx == 0.0) */
/*       SSn[ii] = 0.0; */
/*     else */
/*       SSn[ii] = numer/denom; */

/*   } */

/* } */


/* // ------------------------------------------------------------------- */
/* // FUNCTIONS FOR OUTPUT AND INPUT */
/* // ------------------------------------------------------------------- */

/* // write text files for output */
/* void write_text_dynamic_qstls(double *SSn, double *WW, double *psi_re, */
/* 			      double *psi_im, input in){ */

/*   // Static structure factor */
/*   write_text_dsf(SSn, WW, in); */

/*   FILE* fid; */
  
/*   char out_name[100]; */
/*   sprintf(out_name, "psire_rs%.3f_theta%.3f_x%.3f_%s.dat", in.rs, in.Theta, */
/* 	  in.dyn_xtarget, in.theory); */
/*   fid = fopen(out_name, "w"); */
/*   if (fid == NULL) { */
/*     fprintf(stderr, "Error while creating the output file for the dynamic structure factor\n"); */
/*     exit(EXIT_FAILURE); */
/*   } */
/*   for (int ii = 0; ii < in.nW; ii++) */
/*     fprintf(fid, "%.8e %.8e\n", WW[ii], psi_re[ii]); */
  
/*   fclose(fid); */

  
/*   sprintf(out_name, "psiim_rs%.3f_theta%.3f_x%.3f_%s.dat", in.rs, in.Theta, */
/* 	  in.dyn_xtarget, in.theory); */
/*   fid = fopen(out_name, "w"); */
/*   if (fid == NULL) { */
/*     fprintf(stderr, "Error while creating the output file for the dynamic structure factor\n"); */
/*     exit(EXIT_FAILURE); */
/*   } */
/*   for (int ii = 0; ii < in.nW; ii++) */
/*     fprintf(fid, "%.8e %.8e\n", WW[ii], psi_im[ii]); */
  
/*   fclose(fid); */
  
/* } */
