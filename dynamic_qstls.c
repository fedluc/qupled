#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "solvers.h"
#include "utils.h"
#include "dynamic_stls.h"
#include "dynamic_qstls.h"

// -------------------------------------------------------------------
// LOCAL FUNCTIONS
// -------------------------------------------------------------------

// Allocate and free arrays
static void alloc_dynamic_qstls_2Darrays(input in, double **psi_re,
					 double **psi_im);

static void free_dynamic_qstls_2Darrays(double *psi_re,
					double *psi_im);

// Auxiliary density response
static void compute_dynamic_adr(double *psi_re,  double *psi_im,
				double *WW,  double *SS,
				double *xx, input in);

// Auxiliary density response (real part)
static void compute_dynamic_adr_2D(double *psi_re, double *psi_im,
				   double *WW, double *SS, double *xx,
				   input in);

static void compute_dynamic_adr_re_lev1(double *psi_re,  double *WW,
					double *SS,  double *xx,
					input in);

static double adr_re_lev1_xW(double ww,  void* pp);

static void compute_dynamic_adr_re_lev2(double *int_lev1,  double WW,
					double xx, double *ww, input in);

static double adr_re_lev2_xwW(double uu,  void* pp);

static void compute_dynamic_adr_re_lev3(double *int_lev2,  double WW,
					double ww,  double *qq,  double *uu,
					input in);

static double adr_re_lev3_xwuW(double qq,  void* pp);

static double adr_re_lev3_xwu0(double qq,  void* pp);

// Auxiliary density response (imaginary part)
static void compute_dynamic_adr_im_lev1(double *psi_re,  double *WW,
					double *SS,  double *xx,
					input in);

static double adr_im_lev1_xW(double ww,  void* pp);

static void compute_dynamic_adr_im_lev2(double *psi_im_lev1,  double WW,
					double xx, double *ww, input in);

static double adr_im_lev2_xwW(double uu,  void* pp);

static double adr_im_lev2_xw0(double uu,  void* pp);

static void compute_dynamic_adr_im_lev3(double *psi_im_lev2, double WW,
					double ww,  double *qq,
					double *uu, input in);

static double adr_im_lev3_xwuW(double qq,  void* pp);

// Dynamic structure factor
static void compute_dsf_qstls(double *SSn, double *phi_re,
			      double *phi_im, double *psi_re,
			      double *psi_im, double *WW,
			      input in);

// Input and output
static void write_text_adr(double *psi_re, double *psi_im,
			   double *WW, input in);


static void write_fixed_dynamic_adr(double *psi_re, double *psi_im,
				    input in);


// -------------------------------------------------------------------
// LOCAL CONSTANTS AND DATA STRUCTURES
// -------------------------------------------------------------------

// Number of data points for the integration of the second part
// auxiliary density response
#define ADR_NU 201

// Parameters for the integrals in the auxiliary density response
struct adr_lev1_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *int_lev1_sp_ptr;
  gsl_interp_accel *int_lev1_acc_ptr;

};

struct adr_lev2_params {

  double ww;
  double xx;
  double Theta;
  double mu;
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
// FUNCTION USED TO COMPUTE THE DYNAMIC PROPERTIES OF QSTLS SCHEME
// -------------------------------------------------------------------

void compute_dynamic_qstls(input in, bool verbose) {

  // Arrays 
  double *WW = NULL; 
  double *phi_re = NULL;
  double *phi_im = NULL;
  double *psi_re = NULL;
  double *psi_im = NULL;
  double *SSn = NULL;
  double *SS = NULL;
  double *xx = NULL;
  
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
  
  // Chemical potential and frequency grid
  init_fixed_dynamic_stls_arrays(&in, WW, verbose);

  // Static structure factor
  if (verbose) printf("Static structure factor (from file): ");
  get_ssf(&SS, &xx, &in);
  if (verbose) printf("Done.\n");

  // Ideal density response
  if (verbose) printf("Normalized ideal Lindhard density calculation: ");
  compute_dynamic_idr(phi_re, phi_im, WW, xx, in);
  if (verbose) printf("Done.\n");
 
  // Auxiliary density response
  if (verbose) printf("Auxiliary density calculation: ");
  fflush(stdout);
  compute_dynamic_adr(psi_re, psi_im, WW, SS, xx, in);
  if (verbose) printf("Done.\n");

  // Dynamic structure factor
  if (verbose) printf("Dynamic structure factor calculation: ");
  compute_dsf_qstls(SSn, phi_re, phi_im, psi_re, psi_im, WW, in);
  if (verbose) printf("Done.\n");
  
  // Output to file
  if (verbose) printf("Writing output files: ");
  write_text_dynamic_stls(SSn, phi_re, phi_im, WW, in);
  write_text_dynamic_qstls(psi_re, psi_im, WW, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_dynamic_stls_arrays(WW, phi_re, phi_im, SSn, NULL); //NULL should be replaced with xx
  free_dynamic_qstls_arrays(psi_re, psi_im, SS, xx);
  
 
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_qstls_arrays(input in, double **psi_re, 
			       double **psi_im){

  *psi_re = malloc( sizeof(double) * in.nW);
  if (*psi_re == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the auxiliary density response\n");
    exit(EXIT_FAILURE);
  }
  
  *psi_im = malloc( sizeof(double) * in.nW);
  if (*psi_im == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the auxiliary density response\n");
    exit(EXIT_FAILURE);
  }
  
}

void alloc_dynamic_qstls_2Darrays(input in, double **psi_re,
				 double **psi_im){
 
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

void free_dynamic_qstls_arrays(double *psi_re, double *psi_im,
			       double *SS, double *xx){

  free(psi_re);
  free(psi_im);
  free(SS);
  free(xx);
  
}


void free_dynamic_qstls_2Darrays(double *psi_re, double *psi_im){

  free(psi_re);
  free(psi_im);

}


// ---------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE STATIC STRUCTURE FACTOR (FROM FILE)
// ---------------------------------------------------------------------

void get_ssf(double **SS, double **xx, input *in){

   // Variables
  size_t ssf_file_name_len = 1000;
  char *ssf_file_name;
  input in_tmp = *in;
   
  // File with static structure factor
  if (strcmp(in->dyn_struct_file, NO_FILE_STR)==0){
    ssf_file_name = malloc( sizeof(char) * ssf_file_name_len);
    sprintf(ssf_file_name, "ssf_rs%.3f_theta%.3f_%s.dat",
	    in->rs, in->Theta, in->theory);
  }
  else {
    ssf_file_name_len = strlen(in->dyn_struct_file) + 1;
    ssf_file_name = malloc( sizeof(char) * ssf_file_name_len);
    strcpy(ssf_file_name, in->dyn_struct_file);
  }
 
  // Get size of data stored in the input file
  get_data_format_from_text(ssf_file_name, &in_tmp.nx, &in_tmp.nl);

  // Allocate temporary arrays to store the structural properties
  *SS = malloc( sizeof(double) * in_tmp.nx);
  *xx = malloc( sizeof(double) * in_tmp.nx);
  if (*SS == NULL ||
      *xx == NULL) {
    fprintf(stderr, "Failed to allocate memory for the data read"
  	    " from file\n");
    exit(EXIT_FAILURE);
  }

  // Get data from input file
  get_data_from_text(ssf_file_name, in_tmp.nx, in_tmp.nl,
		     *SS, *xx, &in_tmp);
  in->nx=in_tmp.nx;
  in->dx=in_tmp.dx; // Set by get_restart_data
  in->xmax=in_tmp.xmax; // Set by get_restart_data

  // Free memory
  free(ssf_file_name);
  
}

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE AUXILIARY DENSITY RESPONSE
// ------------------------------------------------------------------

// Auxiliary density response
void compute_dynamic_adr(double *psi_re, double *psi_im,
			 double *WW, double *SS,
			 double *xx, input in) {

  // We compute the auxiliary density response for all the
  // wave-vectors in xx (loaded from the static strucutre
  // factor file). Then we store the all these
  // auxiliary density responses to file for later use and
  // we interpolate them to the wave-vector xTarget
  // in order to proceed with the calculations.

  // Target wave vector
  double xTarget = in.dyn_xtarget;

  // Temporary arrays to store results for multiple wave vectors
  double *psi_re_2D = NULL;
  double *psi_im_2D = NULL;
  double *psi_re_1D = NULL;
  double *psi_im_1D = NULL;

  // Temporary input structure
  input in_1D = in;
  in_1D.nx = 1;
  
  // Variables for interpolation
  gsl_spline *psi_re_sp_ptr;
  gsl_interp_accel *psi_re_acc_ptr;
  gsl_spline *psi_im_sp_ptr;
  gsl_interp_accel *psi_im_acc_ptr;

  // Allocate arrays
  alloc_dynamic_qstls_2Darrays(in, &psi_re_2D, &psi_im_2D);
  alloc_dynamic_qstls_2Darrays(in_1D, &psi_re_1D, &psi_im_1D);
  psi_re_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  psi_re_acc_ptr = gsl_interp_accel_alloc();
  psi_im_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  psi_im_acc_ptr = gsl_interp_accel_alloc();

  // Auxiliary density response for multiple wave vectors
  compute_dynamic_adr_2D(psi_re_2D, psi_im_2D, WW, SS, xx, in);

  // Write the auxiliary density response to file
  write_fixed_dynamic_adr(psi_re_2D, psi_im_2D, in);
  
  // Interpolate to wave-vector given in input
  for (int jj=0; jj<in.nW; jj++){
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
  free_dynamic_qstls_2Darrays(psi_re_2D, psi_im_2D);
  free_dynamic_qstls_2Darrays(psi_re_1D, psi_im_1D);
  gsl_spline_free(psi_re_sp_ptr);
  gsl_interp_accel_free(psi_re_acc_ptr);
  gsl_spline_free(psi_im_sp_ptr);
  gsl_interp_accel_free(psi_im_acc_ptr);

  
}

// Ideal density response (for multiple wave vectors)
void compute_dynamic_adr_2D(double *psi_re, double *psi_im,
			    double *WW, double *SS, double *xx,
			    input in) {

  
  // Real component
  compute_dynamic_adr_re_lev1(psi_re, WW, SS, xx, in);
  
  // Imaginary component
  compute_dynamic_adr_im_lev1(psi_im, WW, SS, xx, in);
  
}

// ------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE REAL PART OF THE AUXILIARY 
// DENSITY RESPONSE
// ------------------------------------------------------------------

// Real part of the auxiliary density response (level 1)
void compute_dynamic_adr_re_lev1(double *psi_re, double *WW,
				  double *SS, double *xx,
				  input in) {

  // Parallel calculations
  #pragma omp parallel
  {
  
    double err;
    size_t nevals;
    double *int_lev1  = malloc( sizeof(double) * in.nx);
    if (int_lev1 == NULL){
      fprintf(stderr, "Failed to allocate memory for calculation"
	      " of the real part of the auxiliary density"
	      " response function\n");
      exit(EXIT_FAILURE);
    }
    
    // Declare accelerator and spline objects
    gsl_spline *ssf_sp_ptr;
    gsl_interp_accel *ssf_acc_ptr;
    gsl_spline *int_lev1_sp_ptr;
    gsl_interp_accel *int_lev1_acc_ptr;
    
    // Allocate the accelerator and the spline objects
    ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    ssf_acc_ptr = gsl_interp_accel_alloc();
    int_lev1_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    int_lev1_acc_ptr = gsl_interp_accel_alloc();
    
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);

    // Interpolate static structure factor
    gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);
    
    // Loop over the frequency
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){
      for (int jj=0; jj<in.nW; jj++){

	// Integration function
	gsl_function ff_int_lev1;
	ff_int_lev1.function = &adr_re_lev1_xW;
	
	// Inner integrals
	compute_dynamic_adr_re_lev2(int_lev1, WW[jj], xx[ii], xx, in);
	gsl_spline_init(int_lev1_sp_ptr, xx, int_lev1, in.nx);
	
	// Integral over w
	struct adr_lev1_params plev1 = {ssf_sp_ptr,
					ssf_acc_ptr,
					int_lev1_sp_ptr,
					int_lev1_acc_ptr};
	ff_int_lev1.params = &plev1;
	gsl_integration_cquad(&ff_int_lev1,
			      xx[0], xx[in.nx-1],
			      0.0, QUAD_REL_ERR,
			      wsp,
			      &psi_re[idx2(ii,jj,in.nx)],
			      &err, &nevals);
      
      }
    }
    
    // Free memory
    free(int_lev1);
    gsl_integration_cquad_workspace_free(wsp);
    gsl_spline_free(ssf_sp_ptr);
    gsl_interp_accel_free(ssf_acc_ptr);
    gsl_spline_free(int_lev1_sp_ptr);
    gsl_interp_accel_free(int_lev1_acc_ptr);
    
  }
  
}

// Integrand for level 1 of the real auxiliary density response (vector = x, frequency = W)
double adr_re_lev1_xW(double ww, void* pp) {
  
  struct adr_lev1_params* params = (struct adr_lev1_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  gsl_spline* int_lev1_sp_ptr = (params->int_lev1_sp_ptr);
  gsl_interp_accel* int_lev1_acc_ptr = (params->int_lev1_acc_ptr);
  double ffp1 = gsl_spline_eval(int_lev1_sp_ptr, ww, int_lev1_acc_ptr);
  double ssfm1 = gsl_spline_eval(ssf_sp_ptr, ww, ssf_acc_ptr) - 1.0;

  return ww*ssfm1*ffp1;

}


// Real part of the auxiliary density response (level 2)
void compute_dynamic_adr_re_lev2(double *int_lev1, double WW,
				 double xx,  double *ww,
				 input in) {

  double err;
  size_t nevals;

  // Integration limits
  double uu[ADR_NU];
  double du = 2.0/(ADR_NU - 1);
  double int_lev2[ADR_NU];
  
  // Declare accelerator and spline objects
  gsl_spline *int_lev2_sp_ptr;
  gsl_interp_accel *int_lev2_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  int_lev2_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, ADR_NU);
  int_lev2_acc_ptr = gsl_interp_accel_alloc();
  
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
     
  // Integration function
  gsl_function ff_int_lev2;
  ff_int_lev2.function = &adr_re_lev2_xwW;
  

  // Fill array with integration variable (u)
  for (int ii=0; ii<ADR_NU; ii++){
    uu[ii] = -1 + du*ii;
  }
  
  // Loop over w (wave-vector)
  for (int ii=0; ii<in.nx; ii++) {

    // Inner integral
    compute_dynamic_adr_re_lev3(int_lev2, WW, ww[ii], ww, uu, in);
    gsl_spline_init(int_lev2_sp_ptr, uu, int_lev2, ADR_NU);
    
    // Integration over u (wave-vector squared)
    struct adr_lev2_params plev2 = {ww[ii], xx,
				    in.Theta, in.mu,
				    int_lev2_sp_ptr,
				    int_lev2_acc_ptr};
    
    ff_int_lev2.params = &plev2;
    gsl_integration_cquad(&ff_int_lev2,
			  uu[0], uu[ADR_NU-1],
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &int_lev1[ii],
			  &err, &nevals);
  }
 
    

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(int_lev2_sp_ptr);
  gsl_interp_accel_free(int_lev2_acc_ptr);
  
}


// Integrand for level 2 of the real auxiliary density response (vectors = {x,w}, frequency = W)
double adr_re_lev2_xwW(double uu, void* pp) {
  
  struct adr_lev2_params* params = (struct adr_lev2_params*)pp;
  double xx = (params->xx);
  double ww = (params->ww);
  gsl_spline* int_lev2_sp_ptr = (params->int_lev2_sp_ptr);
  gsl_interp_accel* int_lev2_acc_ptr = (params->int_lev2_acc_ptr);  
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double denom = xx2 +  ww2 - 2.0*xx*ww*uu;
  double ffp2 = gsl_spline_eval(int_lev2_sp_ptr, uu, int_lev2_acc_ptr);

  return xx*ww*ffp2/denom;
  
}

// Real part of the auxiliary density response (level 3)
void compute_dynamic_adr_re_lev3(double *int_lev2, double WW,
				  double ww, double *qq, double *uu,
				  input in) {

  double err;
  size_t nevals;
  double xx = in.dyn_xtarget;
 
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
    
  // Integration function
  gsl_function ff_int_lev3;
  if (WW == 0.0)
    ff_int_lev3.function = &adr_re_lev3_xwu0;
  else
    ff_int_lev3.function = &adr_re_lev3_xwuW;
  
  // Loop over u (wave-vector squared)
  for (int ii=0; ii<ADR_NU; ii++){
        
    // Integrate over q (wave-vector)
    struct adr_lev3_params plev3 = {in.mu,in.Theta, xx, ww, uu[ii], WW};
    ff_int_lev3.params = &plev3;
    gsl_integration_cquad(&ff_int_lev3,
			  qq[0], qq[in.nx-1],
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
  double qq2 = qq*qq;
  double WW2 = WW*WW;
  double txq = 2.0*xx*qq;
  double tt = xx2 - xx*ww*uu;
  double txqpt = txq + tt;
  double txqmt = txq - tt;
  double txqpt2 = txqpt*txqpt;
  double txqmt2 = txqmt*txqmt;
  double logarg = (txqpt2 - WW2)/(txqmt2 - WW2);

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
  double qq2 = qq*qq;
  double txq = 2.0*xx*qq;
  double tt = xx2 - xx*ww*uu;
  double tt2 = tt*tt;
  double logarg = (tt + txq)/(tt - txq);

  if (xx == 0 || qq == 0){
    return 0;
  }
  else if  (tt == txq){
    return -(3.0/(2.0*Theta))
    *qq2*qq/(exp(qq2/Theta - mu)+ exp(-qq2/Theta + mu) + 2.0);
  }
  else {
    
    if (logarg < 0.0) logarg = -logarg;
    return  -(3.0/(4.0*Theta))
      *qq/(exp(qq2/Theta - mu)+ exp(-qq2/Theta + mu) + 2.0)
      *((qq2 - tt2/(4.0*xx2))*log(logarg) + qq*tt/xx);
				     
  }
  

}


// ------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE IMAGINARY PART OF THE AUXILIARY 
// DENSITY RESPONSE
// ------------------------------------------------------------------

// Imaginary part of the auxiliary density response (level 1)
void compute_dynamic_adr_im_lev1(double *psi_im, double *WW,
				  double *SS, double *xx,
				  input in) {

  // Parallel calculations
  #pragma omp parallel
  {
  
    double err;
    size_t nevals;
    double *int_lev1  = malloc( sizeof(double) * in.nx);
    if (int_lev1 == NULL){
      fprintf(stderr, "Failed to allocate memory for calculation"
	      " of the imaginary part of the auxiliary density"
	      " response function\n");
      exit(EXIT_FAILURE);
    }
    
    // Declare accelerator and spline objects
    gsl_spline *ssf_sp_ptr;
    gsl_interp_accel *ssf_acc_ptr;
    gsl_spline *int_lev1_sp_ptr;
    gsl_interp_accel *int_lev1_acc_ptr;
    
    // Allocate the accelerator and the spline objects
    ssf_sp_ptr = gsl_spline_alloc(gsl_interp_linear, in.nx);
    ssf_acc_ptr = gsl_interp_accel_alloc();
    int_lev1_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    int_lev1_acc_ptr = gsl_interp_accel_alloc();
    
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);
    
    // Loop over the frequency
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){
      for (int jj=0; jj<in.nW; jj++){
	
	// Integration function
	gsl_function ff_int_lev1;
	ff_int_lev1.function = &adr_im_lev1_xW;
	
	// Inner integrals
	compute_dynamic_adr_im_lev2(int_lev1, WW[jj], xx[ii], xx, in);
	
	// Construct integrand
	gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);
	gsl_spline_init(int_lev1_sp_ptr, xx, int_lev1, in.nx);
	
	// Integral over w
	struct adr_lev1_params plev1 = {ssf_sp_ptr,
					ssf_acc_ptr,
					int_lev1_sp_ptr,
					int_lev1_acc_ptr};
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
    free(int_lev1);
    gsl_integration_cquad_workspace_free(wsp);
    gsl_spline_free(int_lev1_sp_ptr);
    gsl_interp_accel_free(int_lev1_acc_ptr);
    
  }
  
}

// Integrand for level 1 of the imaginary auxiliary density response (vector = x, frequency = W)
double adr_im_lev1_xW(double ww, void* pp) {
  
  struct adr_lev1_params* params = (struct adr_lev1_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  gsl_spline* int_lev1_sp_ptr = (params->int_lev1_sp_ptr);
  gsl_interp_accel* int_lev1_acc_ptr = (params->int_lev1_acc_ptr);
  double ffp1 = gsl_spline_eval(int_lev1_sp_ptr, ww, int_lev1_acc_ptr);
  double ssfm1 = gsl_spline_eval(ssf_sp_ptr, ww, ssf_acc_ptr) - 1.0;

  return ww*ssfm1*ffp1;

}


// Imaginary part of the auxiliary density response (level 2)
void compute_dynamic_adr_im_lev2(double *int_lev1, double WW,
				 double xx, double *ww, input in) {

  double err;
  size_t nevals;

  // Integration limits
  double uu[ADR_NU];
  double du = 2.0/(ADR_NU - 1);
  double int_lev2[ADR_NU];
  
  // Declare accelerator and spline objects
  gsl_spline *int_lev2_sp_ptr;
  gsl_interp_accel *int_lev2_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  int_lev2_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, ADR_NU);
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

  // Fill array with integration variable (u)
  for (int ii=0; ii<ADR_NU; ii++){
    uu[ii] = -1 + du*ii;
  }
  
  // Loop over w (wave-vector)
  for (int ii=0; ii<in.nx; ii++) {

    if (WW > 0.0) {

      // Inner integral
      compute_dynamic_adr_im_lev3(int_lev2, WW, ww[ii], ww, uu, in);
      
      // Construct integrand
      gsl_spline_init(int_lev2_sp_ptr, uu, int_lev2, ADR_NU);

    }
    
    // Integration over u (wave-vector squared)
    struct adr_lev2_params plev2 = {ww[ii], xx,
					 in.Theta, in.mu,
					 int_lev2_sp_ptr,
                                         int_lev2_acc_ptr};
    
    ff_int_lev2.params = &plev2;
    gsl_integration_cquad(&ff_int_lev2,
			  uu[0], uu[ADR_NU-1],
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &int_lev1[ii],
			  &err, &nevals);
    
  }
 
    

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(int_lev2_sp_ptr);
  gsl_interp_accel_free(int_lev2_acc_ptr);
  
}


// Integrand for level 2 of the imaginary auxiliary density response (vectors = {x,w}, frequency = W)
double adr_im_lev2_xwW(double uu, void* pp) {
  
  struct adr_lev2_params* params = (struct adr_lev2_params*)pp;
  double xx = (params->xx);
  double ww = (params->ww);
  gsl_spline* int_lev2_sp_ptr = (params->int_lev2_sp_ptr);
  gsl_interp_accel* int_lev2_acc_ptr = (params->int_lev2_acc_ptr);  
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double denom = xx2 +  ww2 - 2.0*xx*ww*uu;
  double ffp2 = gsl_spline_eval(int_lev2_sp_ptr, uu, int_lev2_acc_ptr);

  return (3.0*M_PI/8)*xx*ww*ffp2/denom;
  
}

// Integrand for level 2 of the imaginary auxiliary density response (vectors = {x,w}, frequency = 0)
double adr_im_lev2_xw0(double uu, void* pp) {
  
  struct adr_lev2_params* params = (struct adr_lev2_params*)pp;
  double xx = (params->xx);
  double ww = (params->ww);
  double Theta = (params->Theta);
  double mu = (params->mu);
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double tt = xx2 - xx*ww*uu;
  double tt2 = tt*tt;
  double denom = 2*tt + ww2 - xx2;

  return tt*xx*ww/denom*1.0/(exp(tt2/(4.0*Theta*xx2) - mu) + 1.0);
  
}


// Imaginary part of the auxiliary density response (level 3)
void compute_dynamic_adr_im_lev3(double *int_lev2, double WW,
				 double ww, double *qq, double *uu,
				 input in) {

  double err;
  size_t nevals;
  double q_min;
  double q_max;
  double xx = in.dyn_xtarget;
  double xx2 = xx*xx;
  double xw = xx*ww;
  double x2mxwu;
 
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
    
  // Integration function
  gsl_function ff_int_lev3;
  ff_int_lev3.function = &adr_im_lev3_xwuW;
  
  // Loop over u (wave-vector squared)
  for (int ii=0; ii<ADR_NU; ii++){
    
    // Integration limits
    x2mxwu = xx2 - xw*uu[ii];
    if (x2mxwu < 0.0) x2mxwu = -x2mxwu;
    q_min = (WW - x2mxwu)/(2.0*xx);
    if (q_min < 0.0) q_min = -q_min;
    q_max = (WW + x2mxwu)/(2.0*xx);
    
    // Integrate over q (wave-vector)
    struct adr_lev3_params plev3 = {in.mu,in.Theta, xx, ww, uu[ii], WW};
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

// Integrand for level 3 of the imaginary auxiliary density response (vectors = {x,w,u}, frequency = W)
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
  double hh1 = (xx2 - xx*ww*uu + WW)/(2.0*xx);
  double hh2 = (xx2 - xx*ww*uu - WW)/(2.0*xx);
  double hh12 = hh1*hh1;
  double hh22 = hh2*hh2;
  int out1 = 0;
  int out2 = 0;
  
  if (qq2 > hh12) 
    out1 = 1;

  if (qq2 > hh22)
    out2 = -1;
  
  return (out1 + out2)*qq/(exp(qq2/Theta - mu) + 1.0);

}


// ---------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC STRUCTURE FACTOR
// ---------------------------------------------------------------------

void compute_dsf_qstls(double *SSn, double *phi_re, double *phi_im,
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
  
  for (int ii=0; ii<in.nW; ii++){

    if (WW[ii] == 0.0) {

      ff2 = in.Theta/(4.0*xx);
      numer = (1.0 - ff1*psi_re[ii])
	/(exp(xx*xx/(4.0*in.Theta) - in.mu) + 1)
	- 3.0/(4.0*xx)*ff1*phi_re[ii]*psi_im[ii];
      numer *= ff2;
      denom_re = 1.0 + ff1 * (phi_re[ii] - psi_re[ii]);
      denom = denom_re * denom_re;
	
    }
    else {
      
      ff2 = 1.0/(1.0 - exp(-WW[ii]/in.Theta));
      numer = phi_im[ii] + ff1*(phi_re[ii]*psi_im[ii] -
				phi_im[ii]*psi_re[ii]);
      numer *= (ff2/M_PI);
      denom_re = 1.0 + ff1 * (phi_re[ii] - psi_re[ii]);
      denom_im = ff1 * (phi_im[ii] - psi_im[ii]);
      denom = denom_re*denom_re + denom_im*denom_im;
      
    }
    
    if (xx == 0.0)
      SSn[ii] = 0.0;
    else
      SSn[ii] = numer/denom;

  }

}


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

// write text files for output
void write_text_dynamic_qstls(double *psi_re, double *psi_im,
			      double *WW, input in){

  // Auxiliary density response
  write_text_adr(psi_re, psi_im, WW, in);

}

// write ideal density response
void write_text_adr(double *psi_re, double *psi_im, double *WW, input in){

  FILE* fid;

  char out_name[100];
  sprintf(out_name, "adr_rs%.3f_theta%.3f_x%.3f_%s.dat", in.rs, in.Theta,
	  in.dyn_xtarget, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  for (int ii = 0; ii < in.nW; ii++)
    fprintf(fid, "%.8e %.8e %.8e\n", WW[ii], psi_re[ii], psi_im[ii]);
  
  fclose(fid);
  
}

// write auxiliary density response to binary file
void write_fixed_dynamic_adr(double *psi_re, double *psi_im,
			     input in){
  
  // Name of output file
  char out_name[100];
  sprintf(out_name, "dynamic_restart_rs%.3f_theta%.3f_%s.bin", in.rs, in.Theta,
	  in.theory);

  // Open binary file (we append to the file created for the ideal density response)
  FILE *fid = NULL;
  fid = fopen(out_name, "ab");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file to store the density"
	    " responses\n");
    exit(EXIT_FAILURE);
  }
  
  // Density response
  fwrite(psi_re, sizeof(double), in.nx * in.nW, fid);
  fwrite(psi_im, sizeof(double), in.nx * in.nW, fid);

  // Close binary file
  fclose(fid);

}
