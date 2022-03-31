#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "solvers.h"
#include "utils.h"
#include "stls_iet.h"
#include "dynamic_stls.h"
#include "dynamic_qstls.h"

// -------------------------------------------------------------------
// LOCAL FUNCTIONS
// -------------------------------------------------------------------

// Bridge function term
static void get_bf(double **bf, double *xx, input in);

// Auxiliary density response 
static void compute_dynamic_adr(double *psi_re, double *psi_im,
				double *phi_re, double *phi_im,
				double *WW, double *SS,
				double *bf, double *xx,
				input in);

// Auxiliary density response (partially dynamic) 
static void compute_dynamic_adr_pd(double *phi_re, double *phi_im,
				   double *psi_re, double *psi_im,
				   double *WW, double *SS,
				   double *bf, double *xx,
				   input in);

// Auxiliary density response (fully dynamic) 
static void compute_dynamic_adr_fd(double *phi_re, double *phi_im,
				   double *psi_re, double *psi_im,
				   double *WW, double *SS,
				   double *bf, double *xx,
				   input in);

// Iterations to determine the auxiliary density response
static double adr_err(double *psi_re, double *psi_re_new, input in);

static void adr_update(double *psi_re, double *psi_re_new, input in);

// Auxiliary density response (partially dynamic, real part) 
static void compute_dynamic_adr_re_lev1(double *psi_re_new, double *psi_re,
					double *psi_fixed, double *phi_re,
					double *WW, double *SS, double *bf,
					double *xx, bool compute_fixed,
					input in);

static void compute_dynamic_adr_re_lev1_1(double *int_lev1_1, double *psi_re,
					  double *phi_re, double *SS,
					  double *bf, input in);

static void read_dynamic_adr_fixed(double *var, double *psi_fixed,
				   int ii, int jj, input in);

static void write_dynamic_adr_fixed(double *var, double *psi_fixed,
				    int ii, int jj, input in);

static double adr_re_lev1_xW(double ww, void* pp);

static void compute_dynamic_adr_re_lev2(double *int_lev1, double WW,
					double xx, double *SS,
					double *ww, input in);

static double adr_re_lev2_xwW(double uu, void* pp);

static void compute_dynamic_adr_re_lev3(double *int_lev2, double WW,
					double xx, double ww,
					double *uu, input in);

static double adr_re_lev3_xwuW(double qq, void* pp);

static double adr_re_lev3_xwu0(double qq, void* pp);

// Auxiliary density response (partially dynamic, imaginary part) 
static void compute_dynamic_adr_im_lev1(double *psi_im, double *psi_re,
					double *phi_re, double *WW,
					double *SS, double *bf,
					double *xx, input in);

static void compute_dynamic_adr_im_lev1_1(double *int_lev1_1, double *psi_re,
					  double *phi_re, double *SS,
					  double *bf, input in);

static double adr_im_lev1_xW(double ww, void* pp);

static void compute_dynamic_adr_im_lev2(double *int_lev1, double WW,
					double xx, double *SS,
					double *ww, input in);

static double adr_im_lev2_xwW(double uu, void* pp);

static double adr_im_lev2_xw0(double uu, void* pp);

static void compute_dynamic_adr_im_lev3(double *int_lev2, double WW,
					double xx, double ww,
					double *uu, input in);

static double adr_im_lev3_xwuW(double qq, void* pp);

// Auxiliary density response (fully dynamic, real part) 
static void compute_dynamic_adr_fd_re(double *psi_re_new,
				      double *psi_re, double *psi_im,
				      double *phi_re, double *phi_im,
				      double *psi_fixed_p1,
				      double *psi_fixed_p2,
				      double *WW, double *SS,
				      double *bf, double *xx,
				      bool compute_fixed, input in);

static void compute_dynamic_adr_fd_re_lev1(double *psi_re_new,
					   double *psi_re, double *psi_im,
					   double *phi_re, double *phi_im,
					   double *psi_fixed,
					   double *WW, double *SS,
					   double *bf, double *xx,
					   int ncomp, bool compute_fixed,
					   input in);

static void compute_dynamic_adr_fd_re_lev1_p1_1(double *int_lev1_1, double *psi_re,
						double *psi_im, double *phi_re,
						double *phi_im, double *SS,
						double *bf, int jj, input in);

static void compute_dynamic_adr_fd_re_lev1_p2_1(double *int_lev1_1, double *psi_re,
						double *psi_im, double *phi_re,
						double *phi_im, double *SS,
						int jj, input in);

// Auxiliary density response (fully dynamic, imaginary part) 
static void compute_dynamic_adr_fd_im(double *psi_im_new,
				      double *psi_re, double *psi_im,
				      double *phi_re, double *phi_im,
				      double *psi_fixed_p1,
				      double *psi_fixed_p2,
				      double *psi_fixed_p3,
				      double *WW, double *SS,
				      double *bf, double *xx,
				      bool compute_fixed, input in);

static void compute_dynamic_adr_fd_im_lev1(double *psi_im_new,
					   double *psi_re, double *psi_im,
					   double *phi_re, double *phi_im,
					   double *psi_fixed,
					   double *WW, double *SS,
					   double *bf, double *xx,
					   int ncomp, bool compute_fixed,
					   input in);

static void compute_dynamic_adr_fd_im_zero_frequency_lev1(double *psi_im_new,
							  double *psi_re,
							  double *psi_im,
							  double *phi_re, 
							  double *psi_fixed,
							  double *SS, double *bf,
							  double *xx, int ncomp,
							  bool compute_fixed,
							  input in);

static void compute_dynamic_adr_fd_im_lev1_p1_1(double *int_lev1_1, double *psi_re,
						double *psi_im, double *phi_re,
						double *SS, double *ww,
						double xx, input in);

// Dynamic structure factor
static void compute_dsf_qstls_iet(double *SSn, double *phi_re, double *phi_im,
				  double *psi_re, double *psi_im,
				  double *WW ,input in);

// -------------------------------------------------------------------
// LOCAL CONSTANTS AND DATA STRUCTURES
// -------------------------------------------------------------------

// Parameters for the integrals in auxiliary density response
struct adr_lev1_params {

  gsl_spline *int_lev1_1_sp_ptr;
  gsl_interp_accel *int_lev1_1_acc_ptr;
  gsl_spline *int_lev1_2_sp_ptr;
  gsl_interp_accel *int_lev1_2_acc_ptr;

};


struct adr_lev2_params {

  double mu;
  double Theta;
  double xx;
  double ww;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *int_lev2_sp_ptr;
  gsl_interp_accel *int_lev2_acc_ptr;
  
};


struct adr_lev3_params {

  double mu;
  double Theta;
  double xx;
  double ww;
  double uu;
  double WW;
  
};

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC PROPERTIES OF QSTLS-IET SCHEME
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

  // Allocate arrays
  alloc_dynamic_stls_arrays(in, &WW, &phi_re, &phi_im, &SSn);
  alloc_dynamic_qstls_arrays(in, &psi_re, &psi_im);
  alloc_stls_iet_arrays(in, &bf);
  
  // Chemical potential and frequency grid
  init_fixed_dynamic_stls_arrays(&in, WW, verbose);
 
  // Static structure factor
  if (verbose) printf("Static structure factor (from file): ");
  get_ssf(&SS, &xx, &in);
  if (verbose) printf("Done.\n");

  // Ideal density response
  if (verbose) printf("Ideal density response calculation: ");
  compute_dynamic_idr(phi_re, phi_im, WW, xx, in);
  if (verbose) printf("Done.\n");

  // Bridge function term
  if (verbose) printf("Bridge function Fourier transform: ");
  get_bf(&bf, xx, in);
  if (verbose) printf("Done.\n");
  
  // Auxiliary density response
  if (verbose) printf("Auxiliary density response calculation: ");
  fflush(stdout);
  compute_dynamic_adr(psi_re, psi_im, phi_re, phi_im, WW, SS,
		      bf, xx, in);
  if (verbose) printf("Done.\n");

  // Dynamic structure factor
  if (verbose) printf("Dynamic structure factor calculation: ");
  compute_dsf_qstls_iet(SSn, phi_re, phi_im, psi_re, psi_im, WW, in);
  if (verbose) printf("Done.\n");
  
  // Output to file
  if (verbose) printf("Writing output files: ");
  write_text_dynamic_stls(SSn, phi_re, phi_im, WW, in);
  write_text_dynamic_qstls(psi_re, psi_im, WW, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_dynamic_stls_arrays(WW, phi_re, phi_im, SSn, NULL); // NULL should be replaced with xx
  free_dynamic_qstls_arrays(psi_re, psi_im, SS, xx);
  free_stls_iet_arrays(bf);
  
 
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_qstls_iet_arrays(input in, double **phi_re_2D,
				    double **phi_im_2D,
				    double **psi_re_2D,
				    double **psi_im_2D,
				    double **phi_re_1D,
				    double **phi_im_1D,
				    double **psi_re_1D,
				    double **psi_im_1D){

  *phi_re_2D = malloc( sizeof(double) * in.nx * in.dyn_nW);
  if (*phi_re_2D == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *phi_im_2D = malloc( sizeof(double) * in.nx * in.dyn_nW);
  if (*phi_im_2D == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *psi_re_2D = malloc( sizeof(double) * in.nx * in.dyn_nW);
  if (*psi_re_2D == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the auxiliary density response\n");
    exit(EXIT_FAILURE);
  }
  
  *psi_im_2D = malloc( sizeof(double) * in.nx * in.dyn_nW);
  if (*psi_im_2D == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the auxiliary density response\n");
    exit(EXIT_FAILURE);
  }

  *phi_re_1D = malloc( sizeof(double) * in.nx);
  if (*phi_re_1D == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *phi_im_1D = malloc( sizeof(double) * in.nx);
  if (*phi_im_1D == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *psi_re_1D = malloc( sizeof(double) * in.nx);
  if (*psi_re_1D == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the auxiliary density response\n");
    exit(EXIT_FAILURE);
  }
  
  *psi_im_1D = malloc( sizeof(double) * in.nx);
  if (*psi_im_1D == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the auxiliary density response\n");
    exit(EXIT_FAILURE);
  }
  
}


void free_dynamic_qstls_iet_arrays(double *phi_re_2D,
				   double *phi_im_2D,
				   double *psi_re_2D,
				   double *psi_im_2D,
				   double *phi_re_1D,
				   double *phi_im_1D,
				   double *psi_re_1D,
				   double *psi_im_1D){

  free(phi_re_2D);
  free(phi_im_2D);
  free(psi_re_2D);
  free(psi_im_2D);
  free(phi_re_1D);
  free(phi_im_1D);
  free(psi_re_1D);
  free(psi_im_1D);
  
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
  double *phi_re_tmp = malloc( sizeof(double) * in.dyn_nW);
  double *phi_im_tmp = malloc( sizeof(double) * in.dyn_nW);
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
    compute_dynamic_idr(phi_re_tmp, phi_im_tmp, WW, xx, in);

    // Copy temporary arrays
    for (int jj=0; jj<in.dyn_nW; jj++) {
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

// Auxiliary density response (real and imaginary part)
void compute_dynamic_adr(double *psi_re, double *psi_im,
			 double *phi_re, double *phi_im,
			 double *WW, double *SS,
			 double *bf, double *xx,
			 input in) {

  // For the auxiliary density response we need to know the
  // density response functions on a grid of wave-vectors. This
  // requires to allocate some dedicated 2D arrays and to
  // recompute the ideal density response

  // Target wave vector
  double xTarget = in.dyn_xtarget;

  // Temporary arrays to store results for multiple wave vectors
  double *phi_re_2D = NULL;
  double *phi_im_2D = NULL;
  double *psi_re_2D = NULL;
  double *psi_im_2D = NULL;
  double *phi_re_1D = NULL;
  double *phi_im_1D = NULL;
  double *psi_re_1D = NULL;
  double *psi_im_1D = NULL;

  // Temporary input structure
  input in_1D = in;
  in_1D.dyn_nW = 1;

  // Variables for interpolation
  gsl_spline *psi_re_sp_ptr;
  gsl_interp_accel *psi_re_acc_ptr;
  gsl_spline *psi_im_sp_ptr;
  gsl_interp_accel *psi_im_acc_ptr;
  
  // Allocate arrays
  alloc_dynamic_stls_2Darrays(in, &phi_re_2D, &phi_im_2D);
  alloc_dynamic_qstls_2Darrays(in, &psi_re_2D, &psi_im_2D);
  alloc_dynamic_stls_2Darrays(in_1D, &phi_re_1D, &phi_im_1D);
  alloc_dynamic_qstls_2Darrays(in_1D, &psi_re_1D, &psi_im_1D);
  psi_re_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  psi_re_acc_ptr = gsl_interp_accel_alloc();
  psi_im_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  psi_im_acc_ptr = gsl_interp_accel_alloc();

  // Density response
  if (strcmp(in.dyn_restart_file, NO_FILE_STR)!=0) {
    read_bin_dynamic_adr_2D(psi_re_2D, psi_im_2D, in);      
  }
  else {


    // Ideal
    compute_dynamic_idr_iet(phi_re_2D, phi_im_2D,
			    WW, xx, in);
    
    // Auxiliary (fully or partially dynamic)
    if (in.qstls_iet_static)
      compute_dynamic_adr_pd(phi_re_2D, phi_im_2D,
			     psi_re_2D, psi_im_2D,
			     WW, SS, bf, xx, in);
    else
      compute_dynamic_adr_fd(phi_re_2D, phi_im_2D,
			     psi_re_2D, psi_im_2D,
			     WW, SS, bf, xx, in);

    // Write to file
    write_bin_dynamic_adr_2D(psi_re_2D, psi_im_2D, in);
     
  }

  // Interpolate to wave-vector given in input
  for (int jj=0; jj<in.dyn_nW; jj++){

    for (int ii=0; ii<in.nx; ii++){
      psi_re_1D[ii] = psi_re_2D[idx2(ii,jj,in.nx)];
      psi_im_1D[ii] = psi_im_2D[idx2(ii,jj,in.nx)];
    }
    gsl_spline_init(psi_re_sp_ptr, xx, psi_re_1D, in.nx);
    gsl_spline_init(psi_im_sp_ptr, xx, psi_im_1D, in.nx);
    psi_re[jj] = gsl_spline_eval(psi_re_sp_ptr, xTarget, psi_re_acc_ptr);
    psi_im[jj] = gsl_spline_eval(psi_im_sp_ptr, xTarget, psi_im_acc_ptr);
  }

  // Free memory
  free_dynamic_stls_2Darrays(phi_re_2D, phi_im_2D);
  free_dynamic_qstls_2Darrays(psi_re_2D, psi_im_2D);
  free_dynamic_qstls_2Darrays(phi_re_1D, phi_im_1D);
  free_dynamic_qstls_2Darrays(psi_re_1D, psi_im_1D);
  gsl_spline_free(psi_re_sp_ptr);
  gsl_interp_accel_free(psi_re_acc_ptr);
  gsl_spline_free(psi_im_sp_ptr);
  gsl_interp_accel_free(psi_im_acc_ptr);

}



// Partially dynamic auxiliary density response
void compute_dynamic_adr_pd(double *phi_re, double *phi_im,
			    double *psi_re, double *psi_im,
			    double *WW, double *SS,
			    double *bf, double *xx,
			    input in){

  // Variables for iterative solution
  int iter_counter;
  double iter_err;
  bool compute_fixed;

  // Allocate temporary arrays
  double *psi_re_new = malloc( sizeof(double) * in.nx * in.dyn_nW);
  double *psi_fixed = malloc( sizeof(double) * in.nx * in.dyn_nW * in.nx);
  if (psi_re_new == NULL ||
      psi_fixed == NULL){
    fprintf(stderr, "Failed to allocate memory for calculation"
	    " of the real part of the auxiliary density"
	    " response function\n");
    exit(EXIT_FAILURE);
  }

  // Initialize the auxiliary density response
  for (int ii=0; ii<in.nx; ii++){
    for (int jj=0; jj<in.dyn_nW; jj++){
      psi_re[idx2(ii,jj,in.nx)] = 0.0;
      psi_fixed[idx2(ii,jj,in.nx)] = 0.0;
    }
  }

  // Iterations
  iter_counter = 0;
  iter_err = 1;
  compute_fixed = true;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter) {

    // Update real component
    compute_dynamic_adr_re_lev1(psi_re_new, psi_re, psi_fixed,
				phi_re, WW, SS, bf, xx,
				compute_fixed, in);

    // Compute fixed component only once
    compute_fixed = false;
    
    // Update diagnostics
    iter_counter++;
    iter_err = adr_err(psi_re, psi_re_new, in);
    adr_update(psi_re, psi_re_new, in);

  }

  // Imaginary component
  compute_dynamic_adr_im_lev1(psi_im, psi_re,
			      phi_re, WW, SS,
			      bf, xx, in);
  
  // Free memory
  free(psi_re_new);  
  free(psi_fixed);

}


// Fully dynamic auxiliary density response
void compute_dynamic_adr_fd(double *phi_re, double *phi_im,
			    double *psi_re, double *psi_im,
			    double *WW, double *SS,
			    double *bf, double *xx,
			    input in){

  // Variables for iterative solution
  int iter_counter;
  double iter_err;
  bool compute_fixed;
  
  // Allocate temporary arrays
  double *psi_re_new = malloc( sizeof(double) * in.nx * in.dyn_nW);
  double *psi_im_new = malloc( sizeof(double) * in.nx * in.dyn_nW);
  double *psi_fixed_p1 = malloc( sizeof(double) * in.nx * in.dyn_nW * in.nx);
  double *psi_fixed_p2 = malloc( sizeof(double) * in.nx * in.dyn_nW * in.nx);
  double *psi_fixed_p3 = malloc( sizeof(double) * in.nx * in.dyn_nW * in.nx);
  if (psi_re_new == NULL ||
      psi_im_new == NULL ||
      psi_fixed_p1 == NULL ||
      psi_fixed_p2 == NULL ||
      psi_fixed_p3 == NULL){
    fprintf(stderr, "Failed to allocate memory for calculation"
	    " of the real part of the auxiliary density"
	    " response function\n");
    exit(EXIT_FAILURE);
  }

  // Initialize the auxiliary density response
  for (int ii=0; ii<in.nx; ii++){
    for (int jj=0; jj<in.dyn_nW; jj++){
      psi_re[idx2(ii,jj,in.nx)] = 0.0;
      psi_im[idx2(ii,jj,in.nx)] = 0.0;
      for (int kk=0; kk<in.nx; kk++){
	psi_fixed_p1[idx3(ii,jj,kk,in.nx,in.dyn_nW)] = 0.0;
	psi_fixed_p2[idx3(ii,jj,kk,in.nx,in.dyn_nW)] = 0.0;
	psi_fixed_p3[idx3(ii,jj,kk,in.nx,in.dyn_nW)] = 0.0;
      }
    }
  }

  // Iterations
  iter_counter = 0;
  iter_err = 1;
  compute_fixed = true;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter) {

    // Update real component
    compute_dynamic_adr_fd_re(psi_re_new,
			      psi_re, psi_im,
			      phi_re, phi_im,
			      psi_fixed_p1,
			      psi_fixed_p2,
			      WW, SS, bf,
			      xx, compute_fixed,
			      in); 

    // Update imaginary component
    compute_dynamic_adr_fd_im(psi_im_new,
			      psi_re, psi_im,
			      phi_re, phi_im,
			      psi_fixed_p1,
			      psi_fixed_p2,
			      psi_fixed_p3,
			      WW, SS, bf,
			      xx, compute_fixed,
			      in); 

    // Compute fixed component only once
    compute_fixed = false;
    
    // Update diagnostics
    iter_counter++;
    iter_err = adr_err(psi_im, psi_im_new, in);
    adr_update(psi_re, psi_re_new, in);
    adr_update(psi_im, psi_im_new, in);

  }
  
  // Free memory
  free(psi_re_new);
  free(psi_im_new);
  free(psi_fixed_p1);
  free(psi_fixed_p2);
  free(psi_fixed_p3);
  
}


// Relative error for the iterations
double adr_err(double *psi, double *psi_new, input in){

  double err = 0.0;
  double err_tmp;
  
  for (int ii=0; ii<in.nx; ii++){
    for (int jj=0; jj<in.dyn_nW; jj++){
      err_tmp = psi[idx2(ii,jj,in.nx)] - psi_new[idx2(ii,jj,in.nx)];
      err += err_tmp*err_tmp;
    }
  }
  
  return sqrt(err);
  
}

// Update iterative solution
void adr_update(double *psi, double *psi_new, input in){

  for (int ii=0; ii<in.nx; ii++){
    for (int jj=0; jj<in.dyn_nW; jj++){
      psi[idx2(ii,jj,in.nx)] = in.a_mix * psi_new[idx2(ii,jj,in.nx)]
      + (1 - in.a_mix)*psi[idx2(ii,jj,in.nx)];
    }
  }
}


// ------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE REAL PART OF THE AUXILIARY
// DENSITY RESPONSE (PARTIALLY DYNAMIC)
// ------------------------------------------------------------------


// Real part of the auxiliary density response (level 1)

void compute_dynamic_adr_re_lev1(double *psi_re_new, double *psi_re,
				 double *psi_fixed, double *phi_re,
				 double *WW, double *SS,
				 double *bf, double *xx,
				 bool compute_fixed, input in) {

  // Parallel calculations
  #pragma omp parallel
  {

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
    
    // Loop over the wave-vectors
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){

      // Loop over the frequencies
      for (int jj=0; jj<in.dyn_nW; jj++) {

	// Integration function
	gsl_function ff_int_lev1;
	ff_int_lev1.function = &adr_re_lev1_xW;

	// Integrand (part 1)
	compute_dynamic_adr_re_lev1_1(int_lev1_1, psi_re, phi_re, SS, bf, in);
	gsl_spline_init(int_lev1_1_sp_ptr, xx, int_lev1_1, in.nx);
	  
	// Integrand (part 2)
	if (compute_fixed) {
	  compute_dynamic_adr_re_lev2(int_lev1_2, WW[jj], xx[ii], SS, xx, in);
	  write_dynamic_adr_fixed(int_lev1_2, psi_fixed, ii, jj, in);
	}
	else {
	  read_dynamic_adr_fixed(int_lev1_2, psi_fixed, ii, jj, in);
	}
	gsl_spline_init(int_lev1_2_sp_ptr, xx, int_lev1_2, in.nx);
	      
	// Integral over w
	struct adr_lev1_params plev1 = {int_lev1_1_sp_ptr,
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
void compute_dynamic_adr_re_lev1_1(double *int_lev1_1, double *psi_re,
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

// Read the fixed component of the auxiliary density response integrals 
void read_dynamic_adr_fixed(double *var, double *psi_fixed,
			    int ii, int jj, input in){
  for (int kk=0; kk<in.nx; kk++)
    var[kk] = psi_fixed[idx3(ii,jj,kk,in.nx,in.dyn_nW)];
}


// Write the fixed component of the auxiliary density response integrals 
void write_dynamic_adr_fixed(double *var, double *psi_fixed,
			     int ii, int jj, input in){
  for (int kk=0; kk<in.nx; kk++){
    psi_fixed[idx3(ii,jj,kk,in.nx,in.dyn_nW)] = var[kk];
  }
}


// Integrand for level 1 of the real auxiliary density response (vector = x, frequency = W)
double adr_re_lev1_xW(double ww, void* pp) {
  
  struct adr_lev1_params* params = (struct adr_lev1_params*)pp;
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
void compute_dynamic_adr_re_lev2(double *int_lev1, double WW,
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
  ff_int_lev2.function = &adr_re_lev2_xwW;

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
    compute_dynamic_adr_re_lev3(int_lev2, WW, xx, ww[ii], ww, in);

    // Construct integrand
    gsl_spline_init(int_lev2_sp_ptr, ww, int_lev2, in.nx);
    
    // Integration over u (wave-vector squared)
    struct adr_lev2_params plev2 = {in.mu, in.Theta,
				    xx, ww[ii],
				    ssf_sp_ptr,
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
double adr_re_lev2_xwW(double uu, void* pp) {
  
  struct adr_lev2_params* params = (struct adr_lev2_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  gsl_spline* int_lev2_sp_ptr = (params->int_lev2_sp_ptr);
  gsl_interp_accel* int_lev2_acc_ptr = (params->int_lev2_acc_ptr);
  double ssfm1 = gsl_spline_eval(ssf_sp_ptr, uu, ssf_acc_ptr) - 1.0;
  double fflv2 = gsl_spline_eval(int_lev2_sp_ptr, uu, int_lev2_acc_ptr);

  return uu*ssfm1*fflv2;
  
}

// Real part of the auxiliary density response (level 3)
void compute_dynamic_adr_re_lev3(double *int_lev2, double WW,
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
    ff_int_lev3.function = &adr_re_lev3_xwu0;
  else
    ff_int_lev3.function = &adr_re_lev3_xwuW;
  
  // Loop over u (wave-vector)
  for (int ii=0; ii<in.nx; ii++){

    // Integrate over q (wave-vector)
    struct adr_lev3_params plev3 = {in.mu,in.Theta, xx, ww,
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
double adr_re_lev3_xwuW(double qq, void* pp) {
  
  struct adr_lev3_params* params = (struct adr_lev3_params*)pp;
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
double adr_re_lev3_xwu0(double qq, void* pp) {
  
  struct adr_lev3_params* params = (struct adr_lev3_params*)pp;
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
  double logarg = f1/f2;

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


// ------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE IMAGINARY PART OF THE AUXILIARY
// DENSITY RESPONSE (PARTIALLY DYNAMIC)
// ------------------------------------------------------------------

// Imaginary part of the auxiliary density response (level 1)
void compute_dynamic_adr_im_lev1(double *psi_im, double *psi_re,
				 double *phi_re, double *WW,
				 double *SS, double *bf,
				 double *xx, input in) {

  // Parallel calculations
  #pragma omp parallel
  {

    double err;
    size_t nevals;
    double *int_lev1_1  = malloc( sizeof(double) * in.nx);
    double *int_lev1_2  = malloc( sizeof(double) * in.nx);
    if (int_lev1_1 == NULL ||
	int_lev1_2 == NULL){
      fprintf(stderr, "Failed to allocate memory for calculation"
	      " of the imaginary part of the auxiliary density"
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
    
    // Loop over the wave-vectors
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){

      // Loop over the frequencies
      for (int jj=0; jj<in.dyn_nW; jj++) {

	// Integration function
	gsl_function ff_int_lev1;
	ff_int_lev1.function = &adr_im_lev1_xW;

	// Integrand (part 1)
	compute_dynamic_adr_im_lev1_1(int_lev1_1, psi_re, phi_re, SS, bf, in);
	gsl_spline_init(int_lev1_1_sp_ptr, xx, int_lev1_1, in.nx);
	  
	// Integrand (part 2)
	compute_dynamic_adr_im_lev2(int_lev1_2, WW[jj], xx[ii], SS, xx, in);
	gsl_spline_init(int_lev1_2_sp_ptr, xx, int_lev1_2, in.nx);
	      
	// Integral over w
	struct adr_lev1_params plev1 = {int_lev1_1_sp_ptr,
					       int_lev1_1_acc_ptr,
					       int_lev1_2_sp_ptr,
					       int_lev1_2_acc_ptr};
	ff_int_lev1.params = &plev1;
	gsl_integration_cquad(&ff_int_lev1,
			      xx[0], xx[in.nx-1],
			      0.0, QUAD_REL_ERR,
			      wsp,
			      &psi_im[idx2(ii,jj,in.nx)],
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

// Integrand for level 1 of the imaginary auxiliary density response (part 1)
void compute_dynamic_adr_im_lev1_1(double *int_lev1_1, double *psi_re,
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

// Integrand for level 1 of the imaginary auxiliary density response (vector = x, frequency = W)
double adr_im_lev1_xW(double ww, void* pp) {
  
  struct adr_lev1_params* params = (struct adr_lev1_params*)pp;
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

// Imaginary part of the auxiliary density response (level 2)

void compute_dynamic_adr_im_lev2(double *int_lev1, double WW,
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
	    " of the imaginary part of the auxiliary density"
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
  if (WW == 0.0)
    ff_int_lev2.function = &adr_im_lev2_xw0;
  else
    ff_int_lev2.function = &adr_im_lev2_xwW;

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
    compute_dynamic_adr_im_lev3(int_lev2, WW, xx, ww[ii], ww, in);

    // Construct integrand
    gsl_spline_init(int_lev2_sp_ptr, ww, int_lev2, in.nx);
    
    // Integration over u (wave-vector squared)
    struct adr_lev2_params plev2 = {in.mu, in.Theta,
					   xx, ww[ii],
					   ssf_sp_ptr,
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


// Integrand for level 2 of the imaginary auxiliary density response (vectors = {x,w}, frequency = W)
double adr_im_lev2_xwW(double uu, void* pp) {
  
  struct adr_lev2_params* params = (struct adr_lev2_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  gsl_spline* int_lev2_sp_ptr = (params->int_lev2_sp_ptr);
  gsl_interp_accel* int_lev2_acc_ptr = (params->int_lev2_acc_ptr);
  double ssfm1 = gsl_spline_eval(ssf_sp_ptr, uu, ssf_acc_ptr) - 1.0;
  double fflv2 = gsl_spline_eval(int_lev2_sp_ptr, uu, int_lev2_acc_ptr);

  return uu*ssfm1*fflv2;
  
}


// Integrand for level 2 of the imaginary auxiliary density response (vectors = {x,w}, frequency = 0)
double adr_im_lev2_xw0(double uu, void* pp) {
  
  struct adr_lev2_params* params = (struct adr_lev2_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  double xx = (params->xx);
  double ww = (params->ww);
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  double ssfm1 = gsl_spline_eval(ssf_sp_ptr, uu, ssf_acc_ptr) - 1.0;
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double uu2 = uu*uu;
  double tt = xx2 + ww2 - uu2;

  if (xx == 0.0)
    return 0.0;
  else
    return 0.5*uu*ssfm1*tt/(exp(tt*tt/(16*Theta*xx2) - mu ) + 1);
  
}

// Imaginary part of the auxiliary density response (level 3)

void compute_dynamic_adr_im_lev3(double *int_lev2, double WW,
				     double xx, double ww,
				     double *uu, input in) {

  double err;
  size_t nevals;
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double tt;
  double q_min;
  double q_max;
  
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
    
  // Integration function
  gsl_function ff_int_lev3;
  ff_int_lev3.function = &adr_im_lev3_xwuW;
  
  // Loop over u (wave-vector)
  for (int ii=0; ii<in.nx; ii++){

    // Integration limits
    tt = (xx2 + ww2 - uu[ii]*uu[ii])/2.0;
    if (tt < 0.0) tt = -tt;
    q_min = (WW - tt)/(2.0*xx);
    if (q_min < 0.0) q_min = -q_min;
    q_max = (WW + tt)/(2.0*xx);
    
    
    // Integrate over q (wave-vector)
    struct adr_lev3_params plev3 = {in.mu,in.Theta, xx, ww,
					   uu[ii], WW};
    ff_int_lev3.params = &plev3;
    gsl_integration_cquad(&ff_int_lev3,
			  q_min, q_max,
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &int_lev2[ii],
			  &err, &nevals);
  }

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

// Integrand for level 3 of the imaginary  auxiliary density response (vectors = {x,w,u}, frequency = W)
double adr_im_lev3_xwuW(double qq, void* pp) {
  
  struct adr_lev3_params* params = (struct adr_lev3_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  double xx = (params->xx);
  double ww = (params->ww);
  double uu = (params->uu);
  double WW = (params->WW);
  double qq2 = qq*qq;
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double uu2 = uu*uu;
  double tt = (xx2 + ww2 - uu2)/2.0;
  double hh1 = (tt + WW)/(2.0*xx);
  double hh2 = (tt - WW)/(2.0*xx);
  double hh12 = hh1*hh1;
  double hh22 = hh2*hh2;
  int out1 = 0;
  int out2 = 0;
  
  if (qq2 > hh12)
    out1 = 1;

  if (qq2 > hh22)
    out2 = -1;
  
  return (3.0*M_PI/8.0)*(out1 + out2)*qq/(exp(qq2/Theta - mu) + 1.0);

}



// ------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE REAL PART OF THE AUXILIARY
// DENSITY RESPONSE (FULLY DYNAMIC)
// ------------------------------------------------------------------

// Real part of the auxiliary density response (level 1)
void compute_dynamic_adr_fd_re(double *psi_re_new,
			       double *psi_re, double *psi_im,
			       double *phi_re, double *phi_im,
			       double *psi_fixed_p1,
			       double *psi_fixed_p2,
			       double *WW, double *SS,
			       double *bf, double *xx,
			       bool compute_fixed, input in){

  // Zero frequency contribution
  compute_dynamic_adr_re_lev1(psi_re_new, psi_re,
  			      psi_fixed_p1, phi_re,
  			      WW, SS, bf, xx,
  			      compute_fixed, in);
	 
  // First component (one triple integral)
  compute_dynamic_adr_fd_re_lev1(psi_re_new, psi_re, psi_im,
				 phi_re, phi_im, psi_fixed_p1,
				 WW, SS, bf, xx, 1, compute_fixed,
				 in);

  // Second component (one triple integral)
  compute_dynamic_adr_fd_re_lev1(psi_re_new, psi_re, psi_im,
				 phi_re, phi_im, psi_fixed_p2,
				 WW, SS, bf, xx, 2, compute_fixed,
				 in);
  
  
}

// Real part of the auxiliary density response (one component, level 1)
void compute_dynamic_adr_fd_re_lev1(double *psi_re_new,
				    double *psi_re, double *psi_im,
				    double *phi_re, double *phi_im,
				    double *psi_fixed,
				    double *WW, double *SS,
				    double *bf, double *xx,
				    int ncomp, bool compute_fixed,
				    input in) {

  // Parallel calculations
  #pragma omp parallel
  {

    double err;
    size_t nevals;
    double psi_new_tmp;
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

    // Loop over the wave-vectors
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){

      // Loop over the frequencies (the zero frequency contribution is computed separately)
      for (int jj=1; jj<in.dyn_nW; jj++) {
	
	// Integration function
	gsl_function ff_int_lev1;
	ff_int_lev1.function = &adr_re_lev1_xW;

	// Integrand (part 1)
	switch(ncomp) {
	case 1:
	  compute_dynamic_adr_fd_re_lev1_p1_1(int_lev1_1, psi_re, psi_im,
					      phi_re, phi_im, SS, bf,
					      jj, in);
	  break;
	case 2:
	  compute_dynamic_adr_fd_re_lev1_p2_1(int_lev1_1, psi_re, psi_im,
					      phi_re, phi_im, SS,
					      jj, in);
	  break;
	default:
	  fprintf(stderr, "Unexpected error in the evaluation of the auxiliary"
		  " density response\n");
	  exit(EXIT_FAILURE);
	}
	gsl_spline_init(int_lev1_1_sp_ptr, xx, int_lev1_1, in.nx);
	  
	// Integrand (part 2)
	if (compute_fixed) {
	  
	  switch(ncomp) {
	  case 1:
	    compute_dynamic_adr_re_lev2(int_lev1_2, WW[jj], xx[ii], SS, xx, in);
	    break;
	  case 2:
	    compute_dynamic_adr_im_lev2(int_lev1_2, WW[jj], xx[ii], SS, xx, in);
	    break;
	  }
	  write_dynamic_adr_fixed(int_lev1_2, psi_fixed, ii, jj, in);
	  
	}
	else {
	  read_dynamic_adr_fixed(int_lev1_2, psi_fixed, ii, jj, in);
	}
	gsl_spline_init(int_lev1_2_sp_ptr, xx, int_lev1_2, in.nx);
	      
	// Integral over w
	struct adr_lev1_params plev1 = {int_lev1_1_sp_ptr,
					int_lev1_1_acc_ptr,
					int_lev1_2_sp_ptr,
					int_lev1_2_acc_ptr};
	ff_int_lev1.params = &plev1;
	gsl_integration_cquad(&ff_int_lev1,
			      xx[0], xx[in.nx-1],
			      0.0, QUAD_REL_ERR,
			      wsp,
			      &psi_new_tmp,
			      &err, &nevals);

	// Output
	switch(ncomp) {
	case 1:
	  psi_re_new[idx2(ii,jj,in.nx)] = psi_new_tmp;
	  break;
	case 2:
	  psi_re_new[idx2(ii,jj,in.nx)] -= psi_new_tmp;
	  break;
	}
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

// Integrand for level 1 of the first component of the real auxiliary density response (part 1)
void compute_dynamic_adr_fd_re_lev1_p1_1(double *int_lev1_1, double *psi_re,
					 double *psi_im, double *phi_re,
					 double *phi_im, double *SS,
					 double *bf, int jj, input in){
  
  double numer;
  double denom;
  double phir;
  double phii;
  double psir;
  double psii;
  
  for (int ii=0; ii<in.nx; ii++){
    if (ii==0) int_lev1_1[ii] = 0.0;
    else {
      phir = phi_re[idx2(ii,jj,in.nx)];
      phii = phi_im[idx2(ii,jj,in.nx)];
      psir = psi_re[idx2(ii,jj,in.nx)];
      psii = psi_im[idx2(ii,jj,in.nx)];
      numer = phir*psir + phii*psii;
      denom = phir*phir + phii*phii;
      int_lev1_1[ii] = SS[ii]*(1.0 - bf[ii]) - (numer/denom)*(SS[ii] - 1.0);
    }
  }
  
}


// Integrand for level 1 of the first component of the real auxiliary density response (part 1)
void compute_dynamic_adr_fd_re_lev1_p2_1(double *int_lev1_1, double *psi_re,
					 double *psi_im, double *phi_re,
					 double *phi_im, double *SS,
					 int jj, input in){
  
  double numer;
  double denom;
  double phir;
  double phii;
  double psir;
  double psii;
  
  for (int ii=0; ii<in.nx; ii++){
    if (ii==0) int_lev1_1[ii] = 0.0;
    else {
      phir = phi_re[idx2(ii,jj,in.nx)];
      phii = phi_im[idx2(ii,jj,in.nx)];
      psir = psi_re[idx2(ii,jj,in.nx)];
      psii = psi_im[idx2(ii,jj,in.nx)];
      numer = phii*psir - phir*psii;
      denom = phir*phir + phii*phii;
      int_lev1_1[ii] = (numer/denom)*(SS[ii] - 1.0);
    }
  }
  
}


// ------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE IMAGINARY PART OF THE AUXILIARY
// DENSITY RESPONSE (FULLY DYNAMIC)
// ------------------------------------------------------------------

// Imaginary part of the auxiliary density response (level 1)
void compute_dynamic_adr_fd_im(double *psi_im_new,
			       double *psi_re, double *psi_im,
			       double *phi_re, double *phi_im,
			       double *psi_fixed_p1,
			       double *psi_fixed_p2,
			       double *psi_fixed_p3,
			       double *WW, double *SS,
			       double *bf, double *xx,
			       bool compute_fixed, input in){


  // Zero-frequency contribution, first component
  compute_dynamic_adr_fd_im_zero_frequency_lev1(psi_im_new,
  						psi_re, psi_im,
  						phi_re, psi_fixed_p1,
  						SS, bf, xx, 1,
  						compute_fixed,
  						in);

  // Zero frequency contribution, second component
  compute_dynamic_adr_fd_im_zero_frequency_lev1(psi_im_new,
  						psi_re, psi_im,
  						phi_re, psi_fixed_p3,
  						SS, bf, xx, 2,
  						compute_fixed,
  						in);

  // Finite frequency contribution, first component
  compute_dynamic_adr_fd_im_lev1(psi_im_new, psi_re, psi_im,
				 phi_re, phi_im, psi_fixed_p1,
				 WW, SS, bf, xx, 1, compute_fixed,
				 in);
  
  // Finite frequency contribution, second component
  compute_dynamic_adr_fd_im_lev1(psi_im_new, psi_re, psi_im,
				 phi_re, phi_im, psi_fixed_p2,
				 WW, SS, bf, xx, 2, compute_fixed,
				 in);
  
}


// Imaginary part of the auxiliary density response (one component, level 1)
void compute_dynamic_adr_fd_im_lev1(double *psi_im_new,
				    double *psi_re, double *psi_im,
				    double *phi_re, double *phi_im,
				    double *psi_fixed,
				    double *WW, double *SS,
				    double *bf, double *xx,
				    int ncomp, bool compute_fixed,
				    input in) {

  // Parallel calculations
  #pragma omp parallel
  {

    double err;
    size_t nevals;
    double psi_new_tmp;
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
    
    // Loop over the wave-vectors
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){

      // Loop over the frequencies (the zero frequency contribution is computed separately)
      for (int jj=1; jj<in.dyn_nW; jj++) {
	
	// Integration function
	gsl_function ff_int_lev1;
	ff_int_lev1.function = &adr_re_lev1_xW;

	// Integrand (part 1)
	switch(ncomp) {
	case 1:
	  compute_dynamic_adr_fd_re_lev1_p2_1(int_lev1_1, psi_re, psi_im,
					      phi_re, phi_im, SS,
					      jj, in);
	  break;
	case 2:
	  compute_dynamic_adr_fd_re_lev1_p1_1(int_lev1_1, psi_re, psi_im,
					      phi_re, phi_im, SS, bf,
					      jj, in);
	  break;
	default:
	  fprintf(stderr, "Unexpected error in the evaluation of the auxiliary"
		  " density response\n");
	  exit(EXIT_FAILURE);
	}
	gsl_spline_init(int_lev1_1_sp_ptr, xx, int_lev1_1, in.nx);
	  
	// Integrand (part 2)
	if (compute_fixed){

	  switch(ncomp) {
	  case 1:
	    compute_dynamic_adr_re_lev2(int_lev1_2, WW[jj], xx[ii], SS, xx, in);
	    break;
	  case 2:
	    compute_dynamic_adr_im_lev2(int_lev1_2, WW[jj], xx[ii], SS, xx, in);
	    break;
	  }
	  write_dynamic_adr_fixed(int_lev1_2, psi_fixed, ii, jj, in);
	  
	}
	else {
	  read_dynamic_adr_fixed(int_lev1_2, psi_fixed, ii, jj, in);
	}
	gsl_spline_init(int_lev1_2_sp_ptr, xx, int_lev1_2, in.nx);
	      
	// Integral over w
	struct adr_lev1_params plev1 = {int_lev1_1_sp_ptr,
					int_lev1_1_acc_ptr,
					int_lev1_2_sp_ptr,
					int_lev1_2_acc_ptr};
	ff_int_lev1.params = &plev1;
	gsl_integration_cquad(&ff_int_lev1,
			      xx[0], xx[in.nx-1],
			      0.0, QUAD_REL_ERR,
			      wsp,
			      &psi_new_tmp,
			      &err, &nevals);

	// Output
	switch(ncomp) {
	case 1:
	  psi_im_new[idx2(ii,jj,in.nx)] = psi_new_tmp;
	  break;
	case 2:
	  psi_im_new[idx2(ii,jj,in.nx)] += psi_new_tmp;
	  break;
	}
	
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


// Imaginary part of the auxiliary density response (one component, level 1)
void compute_dynamic_adr_fd_im_zero_frequency_lev1(double *psi_im_new,
						   double *psi_re, double *psi_im,
						   double *phi_re,
						   double *psi_fixed,
						   double *SS, double *bf,
						   double *xx, int ncomp,
						   bool compute_fixed,
						   input in){

  // Parallel calculations
  #pragma omp parallel
  {

    double err;
    size_t nevals;
    double psi_new_tmp;
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
    
    // Loop over the wave-vectors
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){

      // Loop over the frequencies
      for (int jj=0; jj<1; jj++) {
	
	// Integration function
	gsl_function ff_int_lev1;
	ff_int_lev1.function = &adr_re_lev1_xW;

	// Integrand (part 1)
	switch(ncomp) {
	case 1:
	  compute_dynamic_adr_fd_im_lev1_p1_1(int_lev1_1, psi_re,
					      psi_im, phi_re,
					      SS, xx, xx[ii], in);
	  break;
	case 2:
	  compute_dynamic_adr_im_lev1_1(int_lev1_1, psi_re,
					phi_re, SS, bf, in);
	  break;
	default:
	  fprintf(stderr, "Unexpected error in the evaluation of the auxiliary"
		  " density response\n");
	  exit(EXIT_FAILURE);
	}
	gsl_spline_init(int_lev1_1_sp_ptr, xx, int_lev1_1, in.nx);
	  
	// Integrand (part 2)
	if (compute_fixed){

	  switch(ncomp) {
	  case 1:
	    compute_dynamic_adr_re_lev2(int_lev1_2, 0.0, xx[ii], SS, xx, in);
	    break;
	  case 2:
	    compute_dynamic_adr_im_lev2(int_lev1_2, 0.0, xx[ii], SS, xx, in);
	    break;
	  }
	  write_dynamic_adr_fixed(int_lev1_2, psi_fixed, ii, jj, in);
	  
	}
	else
	  read_dynamic_adr_fixed(int_lev1_2, psi_fixed, ii, jj, in);
	gsl_spline_init(int_lev1_2_sp_ptr, xx, int_lev1_2, in.nx);
	      
	// Integral over w
	struct adr_lev1_params plev1 = {int_lev1_1_sp_ptr,
					int_lev1_1_acc_ptr,
					int_lev1_2_sp_ptr,
					int_lev1_2_acc_ptr};
	ff_int_lev1.params = &plev1;
	gsl_integration_cquad(&ff_int_lev1,
			      xx[0], xx[in.nx-1],
			      0.0, QUAD_REL_ERR,
			      wsp,
			      &psi_new_tmp,
			      &err, &nevals);

	// Output
	switch(ncomp) {
	case 1:
	  psi_im_new[idx2(ii,jj,in.nx)] = psi_new_tmp;
	  break;
	case 2:
	  psi_im_new[idx2(ii,jj,in.nx)] += psi_new_tmp;
	  break;
	}
	
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


// Integrand for level 1 of the first component of the imaginary auxiliary density response (part 1)
void compute_dynamic_adr_fd_im_lev1_p1_1(double *int_lev1_1, double *psi_re,
					 double *psi_im, double *phi_re,
					 double *SS, double *ww,
					 double xx, input in){
  
  double ff;
  double xx2 = xx*xx;
  double ww2;
  double phir;
  double psir;
  double psii;
  double phir2;
   
  for (int ii=0; ii<in.nx; ii++){
    if (ii==0) int_lev1_1[ii] = 0.0;
    else {
      ww2 = ww[ii]*ww[ii];
      phir = phi_re[idx2(ii,0,in.nx)];
      psir = psi_re[idx2(ii,0,in.nx)];
      psii = psi_im[idx2(ii,0,in.nx)];
      phir2 = phir*phir;
      ff = ww[ii] * (exp(ww2/(4.0*in.Theta) - in.mu) + 1.0);
      int_lev1_1[ii] = -(4.0/3.0) * xx2 * (SS[ii] - 1.0) * (psir/(phir2*ff) - 4.0/M_PI*psii/phir );
      // In the above -(4.0/3.0) compensates for the -3.0/(8.0) factor in
      // compute_dynamic_adr_re_lev2 and includes also an additional 1/2
      // factor for the correct calculation of the dynamic structure factor
      // (see, compute_dsf_qstls_iet)
    }
  }
  
}

// ---------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC STRUCTURE FACTOR
// ---------------------------------------------------------------------

void compute_dsf_qstls_iet(double *SSn, double *phi_re, double *phi_im,
			   double *psi_re, double *psi_im,
			   double *WW, input in){
    
  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double xx = in.dyn_xtarget;
  double ff1 = 4.0*lambda*in.rs/(M_PI*xx*xx);
  double ff2;
  double numer;
  double denom;
  double denom_re;
  double denom_im;
  double bb[1];
  double xb[1];
  input in_tmp = in;

  // Bridge function at the target wave-vector
  xb[0] = xx;
  in_tmp.nx = 1;
  compute_bridge_function(bb, xb, in_tmp);
 
  // Dynamic structure factor
  // NOTE: The zero frequency contribution of the imaginary
  // part of the auxiliary density response (psi_im[0])
  // is used to store its frequency derivative
  for (int ii=0; ii<in.dyn_nW; ii++){

    if (WW[ii] == 0.0) {

      ff2 = in.Theta/(4.0*xx);
      numer = (1.0 - ff1*psi_re[ii])
    	/(exp(xx*xx/(4.0*in.Theta) - in.mu) + 1)
    	- 3.0/(4.0*xx)*ff1*phi_re[ii]*psi_im[ii];
      numer *= ff2;
      denom_re = 1.0 + ff1 * ((1 - bb[0])*phi_re[ii] - psi_re[ii]);
      denom = denom_re * denom_re;
	
    }
    else {
      
      ff2 = 1.0/(1.0 - exp(-WW[ii]/in.Theta));
      numer = phi_im[ii] + ff1*(phi_re[ii]*psi_im[ii] -
				phi_im[ii]*psi_re[ii]);
      numer *= (ff2/M_PI);
      denom_re = 1.0 + ff1 * ((1 - bb[0])*phi_re[ii] - psi_re[ii]);
      denom_im = ff1 * ((1 - bb[0])*phi_im[ii] - psi_im[ii]);
      denom = denom_re*denom_re + denom_im*denom_im;
      
    }
        
    if (xx == 0.0)
      SSn[ii] = 0.0;
    else
      SSn[ii] = numer/denom;

  }

}
