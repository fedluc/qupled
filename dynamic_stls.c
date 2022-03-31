#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "solvers.h"
#include "utils.h"
#include "chemical_potential.h"
#include "dynamic_stls.h"

// -------------------------------------------------------------------
// LOCAL FUNCTIONS
// -------------------------------------------------------------------

// Frequency grid
static void frequency_grid(double *WW, input *in);

// Ideal density response
static void compute_dynamic_idr_2D(double *phi_re, double *phi_im,
				   double *WW, double *xx,
				   input in);

static void compute_dynamic_idr_re(double *phi_re, double *WW,
				   double *xx, input in);

static void compute_dynamic_idr_im(double *phi_im, double *WW,
				   double *xx, input in);

static double idr_re_xw(double yy, void *pp);

static double idr_re_x0(double yy, void *pp);

static double idr_im_xw(double yy, void *pp);

static double idr_im_x0(double yy, void *pp);

// Static local field correction (from file)
static void get_slfc(double *GG, double **xx, input in);

// Dynamic structure factor
static void compute_dsf(double *SSn, double *phi_re, double *phi_im,
			double GG, double *WW, input in);

// Intermediate scattering function
static void compute_isf(double *FF, double *tt, double *SSn,
			double *WW, input in);
 
static double isf(double WW, void *pp);

// Input and output
static void write_text_dsf(double *SSn, double *WW, input in);

static void write_text_isf(double *SSn, double *ww, input in);

static void write_text_idr(double *phi_re, double *phi_im,
			   double *WW, input in);

static void write_bin_idr_2D(double *phi_re, double *phi_im,
			     input in);

static void read_bin_idr_2D(double *phi_re, double *phi_im,
			    input in);


// -------------------------------------------------------------------
// LOCAL CONSTANTS AND DATA STRUCTURES
// -------------------------------------------------------------------

// Number of data points for imaginary time (varies between 0 and 1)
#define ISF_NTAU 100

// Parameters for integrals in the ideal density response
struct idr_params {

  double xx;
  double mu;
  double Theta;
  double WW;

};

// Parameters for the integrals in the intermediate scattering function
struct isf_params {

  double Theta;
  double tau;
  gsl_spline *dsf_sp_ptr;
  gsl_interp_accel *dsf_acc_ptr;
  
};


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC PROPERTIES OF THE CLASSICAL
// SCHEMES (STLS, VS-STLS AND STLS-IET)
// -------------------------------------------------------------------

void compute_dynamic_stls(input in, bool verbose) {

  // Arrays 
  double *WW = NULL; 
  double *phi_re = NULL;
  double *phi_im = NULL;
  double *SSn = NULL;
  double *xx = NULL;
  
  // Scalars
  double GG;
 
  // Safeguard
  if (in.Theta == 0) {
    printf("Ground state calculations of the dynamic properties"
	   " are not yet implemented.");
    exit(EXIT_FAILURE);
  }
      
  // Get the size of the frequency grid
  get_frequency_grid_size(&in);
  
  // Allocate arrays
  alloc_dynamic_stls_arrays(in, &WW, &phi_re, &phi_im,
			    &SSn);

  // Chemical potential and frequency grid
  init_fixed_dynamic_stls_arrays(&in, WW, verbose);

  /* // Static local field correction (this sets xx) */
  if (verbose) printf("Static local field correction (from file): ");
  get_slfc(&GG, &xx, in);
  if (verbose) printf("Done.\n");
  
  // Ideal density response
  if (verbose) printf("Ideal density response calculation: ");
  compute_dynamic_idr(phi_re, phi_im, WW, xx, in);
  if (verbose) printf("Done.\n");
  
  // Dynamic structure factor
  if (verbose) printf("Dynamic structure factor calculation: ");
  compute_dsf(SSn, phi_re, phi_im, GG, WW, in);
  if (verbose) printf("Done.\n");
  
  /* // Output to file */
  if (verbose) printf("Writing output files: ");
  write_text_dynamic_stls(SSn, phi_re, phi_im, WW, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_dynamic_stls_arrays(WW, phi_re, phi_im, SSn, xx);

 
}

// -------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE SIZE OF THE FREQUENCY GRID
// -------------------------------------------------------------------

void get_frequency_grid_size(input *in){

  // Number of grid points in the frequency grid
  in->dyn_nW = (int)floor(in->dyn_Wmax/in->dyn_dW);
  
}
// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_stls_arrays(input in, double **WW, double **phi_re, 
			       double **phi_im, double **SSn){

  *WW = malloc( sizeof(double) * in.dyn_nW);
  if (*WW == NULL) {
    fprintf(stderr, "Failed to allocate memory for the frequency grid\n");
    exit(EXIT_FAILURE);
  }

  *phi_re = malloc( sizeof(double) * in.dyn_nW);
  if (*phi_re == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *phi_im = malloc( sizeof(double) * in.dyn_nW);
  if (*phi_im == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *SSn = malloc( sizeof(double) * in.dyn_nW);
  if (*SSn == NULL) {
    fprintf(stderr, "Failed to allocate memory for the dynamic structure factor\n");
    exit(EXIT_FAILURE);
  }

}

void alloc_dynamic_stls_2Darrays(input in, double **phi_re,
				 double **phi_im){
 
  *phi_re = malloc( sizeof(double) * in.nx * in.dyn_nW);
  if (*phi_re == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *phi_im = malloc( sizeof(double) * in.nx * in.dyn_nW);
  if (*phi_im == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
 
}


void free_dynamic_stls_arrays(double *WW, double *phi_re, 
			      double *phi_im, double *SSn,
			      double *xx){

  free(WW);
  free(phi_re);
  free(phi_im);
  free(SSn);
  free(xx);
 
}

void free_dynamic_stls_2Darrays(double *phi_re, double *phi_im){

  free(phi_re);
  free(phi_im);

}

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_dynamic_stls_arrays(input *in, double *WW, bool verbose){

  // Print on screen the parameter used to solve the STLS equation
  printf("------ Parameters used in the solution -------------\n");
  printf("Quantum degeneracy parameter: %f\n", in->Theta);
  printf("Quantum coupling parameter: %f\n", in->rs);
  printf("Chemical potential (low and high bound): %f %f\n", 
	 in->mu_lo, in->mu_hi);
  printf("Target wave-vector: %f\n", in->dyn_xtarget);
  printf("Frequency cutoff: %f\n", in->dyn_Wmax);
  printf("Frequency resolution: %f\n", in->dyn_dW);
  printf("----------------------------------------------------\n");
 
  // Chemical potential
  if (in->Theta > 0) {
    if (verbose) printf("Chemical potential calculation: ");
    in->mu = compute_chemical_potential(*in);
    if (verbose) printf("Done. Chemical potential: %.8f\n", in->mu);
  }
  
  // Frequency grid
  if (verbose) printf("Frequency grid initialization: ");
  frequency_grid(WW, in);
  if (verbose) printf("Done.\n");

}

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE FREQUENCY GRID
// ------------------------------------------------------------------

void frequency_grid(double *WW, input *in){
 
  WW[0] = 0.0;
  for (int ii=1; ii < in->dyn_nW; ii++) WW[ii] = WW[ii-1] + in->dyn_dW;

}

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE IDEAL DENSITY RESPONSE
// ------------------------------------------------------------------

// Ideal density response
void compute_dynamic_idr(double *phi_re, double *phi_im,
			 double *WW, double *xx, input in) {


  // We compute the ideal density response for all the
  // wave-vectors in xx (loaded from the static local
  // field correction file). Then we store the all these
  // ideal density responses to file for later use and
  // we interpolate them to the wave-vector xTarget
  // in order to proceed with the calculations.

  // Target wave vector
  double xTarget = in.dyn_xtarget;

  // Temporary arrays to store results for multiple wave vectors
  double *phi_re_2D = NULL;
  double *phi_im_2D = NULL;
  double *phi_re_1D = NULL;
  double *phi_im_1D = NULL;

  // Temporary input structure
  input in_1D = in;
  in_1D.dyn_nW = 1;
  
  // Variables for interpolation
  gsl_spline *phi_re_sp_ptr;
  gsl_interp_accel *phi_re_acc_ptr;
  gsl_spline *phi_im_sp_ptr;
  gsl_interp_accel *phi_im_acc_ptr;

  // Allocate arrays
  alloc_dynamic_stls_2Darrays(in, &phi_re_2D, &phi_im_2D);
  alloc_dynamic_stls_2Darrays(in_1D, &phi_re_1D, &phi_im_1D);
  
  phi_re_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  phi_re_acc_ptr = gsl_interp_accel_alloc();
  phi_im_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  phi_im_acc_ptr = gsl_interp_accel_alloc();

  // Ideal density response for multiple wave vectors
  if (strcmp(in.dyn_restart_file, NO_FILE_STR) != 0) {
    read_bin_idr_2D(phi_re_2D, phi_im_2D, in);
  }
  else{
    compute_dynamic_idr_2D(phi_re_2D, phi_im_2D,
  			   WW, xx, in);
    write_bin_idr_2D(phi_re_2D, phi_im_2D, in);
  }
  
  // Interpolate to wave-vector given in input
  for (int jj=0; jj<in.dyn_nW; jj++){
    for (int ii=0; ii<in.nx; ii++){
      phi_re_1D[ii] = phi_re_2D[idx2(ii,jj,in.nx)];
      phi_im_1D[ii] = phi_im_2D[idx2(ii,jj,in.nx)];
    }
    gsl_spline_init(phi_re_sp_ptr, xx, phi_re_1D, in.nx);
    gsl_spline_init(phi_im_sp_ptr, xx, phi_im_1D, in.nx);
    phi_re[jj] = gsl_spline_eval(phi_re_sp_ptr, xTarget, phi_re_acc_ptr);
    phi_im[jj] = gsl_spline_eval(phi_im_sp_ptr, xTarget, phi_im_acc_ptr);
  }

  // Free memory
  free_dynamic_stls_2Darrays(phi_re_2D, phi_im_2D);
  free_dynamic_stls_2Darrays(phi_re_1D, phi_im_1D);
  gsl_spline_free(phi_re_sp_ptr);
  gsl_interp_accel_free(phi_re_acc_ptr);
  gsl_spline_free(phi_im_sp_ptr);
  gsl_interp_accel_free(phi_im_acc_ptr);

}

// Ideal density response (for multiple wave vectors)
void compute_dynamic_idr_2D(double *phi_re, double *phi_im,
			    double *WW, double *xx, input in) {

  
  // Real component
  compute_dynamic_idr_re(phi_re, WW, xx, in);

  // Imaginary component
  compute_dynamic_idr_im(phi_im, WW, xx, in);
  
}

// Real part of the ideal density response
void compute_dynamic_idr_re(double *phi_re, double *WW,
			    double *xx, input in) {

  double err;
  size_t nevals;
  
  // Integration workspace 
  gsl_integration_cquad_workspace *wsp 
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  if (WW == 0) ff_int.function = &idr_re_x0;
  else ff_int.function = &idr_re_xw;

  // Normalized ideal Lindhard density
  for (int ii=0; ii<in.nx; ii++) {
    for (int jj=0; jj<in.dyn_nW; jj++) {

      struct idr_params phixwp = {xx[ii], in.mu, in.Theta, WW[jj]};
      ff_int.params = &phixwp;
      gsl_integration_cquad(&ff_int, 
			    0.0, in.xmax, 
			    0.0, QUAD_REL_ERR, 
			    wsp, 
			    &phi_re[idx2(ii,jj,in.nx)],
			    &err, &nevals);

    }
  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

// Imaginary part of the ideal density response
void compute_dynamic_idr_im(double *phi_im, double *WW,
			    double *xx, input in) {
 
  double ymin;
  double ymax;
  double err;
  size_t nevals;


  // Integration workspace 
  gsl_integration_cquad_workspace *wsp 
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  if (WW == 0) ff_int.function = &idr_im_x0;
  else ff_int.function = &idr_im_xw;

  // Normalized ideal Lindhard density
  for (int ii=0; ii<in.nx; ii++){
    for (int jj=0; jj<in.dyn_nW; jj++) {
      
      if (xx[ii] == 0.0){
	phi_im[idx2(ii,jj,in.nx)] = 0.0;
	continue;
      }
      
      ymin = (xx[ii]/2.0) - WW[jj]/(2.0*xx[ii]);
      if (ymin < 0.0) ymin = -ymin;
      ymax = (xx[ii]/2.0) + WW[jj]/(2.0*xx[ii]);
      
      struct idr_params phixwp = {xx[ii], in.mu, in.Theta, WW[jj]};
      ff_int.params = &phixwp;
      gsl_integration_cquad(&ff_int, 
			    ymin, ymax, 
			    0.0, QUAD_REL_ERR, 
			    wsp, 
			    &phi_im[idx2(ii,jj,in.nx)],
			    &err, &nevals);
      
    }
  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

// Partial real part of the ideal density response (frequency = w, vector = x)
double idr_re_xw(double yy, void *pp) {

  struct idr_params *params = (struct idr_params*)pp;
  double xx = (params->xx);
  double mu = (params->mu);
  double Theta = (params->Theta);
  double WW = (params->WW);
  double yy2 = yy*yy;
  double ymx = yy - 0.5*xx;
  double ypx = yy + 0.5*xx;
  double ymx2 = ymx*ymx;
  double ypx2 = ypx*ypx;
  double w_2x = WW/(2.0*xx);
  double w_2x2 = w_2x*w_2x;
  double log_arg = (ymx2 - w_2x2)/(ypx2 - w_2x2);


  if (log_arg<0) log_arg = -log_arg;
  
  if (xx > 0.0) {
    return -1.0/(2*xx)*yy/(exp(yy2/Theta - mu) + 1.0)
      *log(log_arg);
  }
  else {
    return 0;
  }

}


// Partial real part of the ideal density response (frequency = 0, vector = x)
double idr_re_x0(double yy, void *pp) {

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


// Partial imaginary part of the ideal density response (frequency = w, vector = x)
double idr_im_xw(double yy, void *pp) {

  struct idr_params *params = (struct idr_params*)pp;
  double xx = (params->xx);
  double mu = (params->mu);
  double Theta = (params->Theta);
  double WW = (params->WW);
  double yy2 = yy*yy;
  double xpw = 0.5*xx + 0.5*WW/xx;
  double xmw = 0.5*xx - 0.5*WW/xx;
  double xpw2 = xpw*xpw;
  double xmw2 = xmw*xmw;
  double out1 = 0;
  double out2 = 0;

  if (xx == 0.0) return 0;

  if (yy2 > xpw2) 
    out1 = 1;

  if (yy2 > xmw2)
    out2 = -1;
  
  return -M_PI/(2*xx)*yy/(exp(yy2/Theta - mu) + 1.0)
    *(out1 + out2);

}

// Partial imaginary part of the ideal density response (frequency = 0, vector = x)
double idr_im_x0(double yy, void *pp) {
  return 0;
}


// ---------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE STATIC LOCAL FIELD CORRECTION (FROM FILE)
// ---------------------------------------------------------------------

void get_slfc(double *GG, double **xx, input in){

  // Variables
  size_t slfc_file_name_len = 1000;
  char *slfc_file_name;
  double *GG_file = NULL;
  gsl_spline *slfc_sp_ptr;
  gsl_interp_accel *slfc_acc_ptr;

  // File with static local field correction
  if (strcmp(in.dyn_struct_file, NO_FILE_STR)==0){
    slfc_file_name = malloc( sizeof(char) * slfc_file_name_len);
    sprintf(slfc_file_name, "slfc_rs%.3f_theta%.3f_%s.dat",
  	    in.rs, in.Theta, in.theory);
  }
  else {
    slfc_file_name_len = strlen(in.dyn_struct_file) + 1;
    slfc_file_name = malloc( sizeof(char) * slfc_file_name_len);
    strcpy(slfc_file_name, in.dyn_struct_file);
  }
 
  // Get size of data stored in the input file
  get_data_format_from_text(slfc_file_name, &in.nx, &in.nl);

  // Allocate temporary arrays to store the structural properties
  GG_file = malloc( sizeof(double) * in.nx);
  *xx = malloc( sizeof(double) * in.nx);
  if (GG_file == NULL ||
      *xx == NULL) {
    fprintf(stderr, "Failed to allocate memory for the data read"
  	    " from file\n");
    exit(EXIT_FAILURE);
  }

  // Get data from input file
  get_data_from_text(slfc_file_name, in.nx, in.nl,
  		     GG_file, *xx, &in);

  // Static local field correction for the wave-vector given in input
  slfc_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  slfc_acc_ptr = gsl_interp_accel_alloc();
  gsl_spline_init(slfc_sp_ptr, *xx, GG_file, in.nx);
  *GG = gsl_spline_eval(slfc_sp_ptr, in.dyn_xtarget, slfc_acc_ptr);
  
  // Free memory
  free(slfc_file_name);
  free(GG_file);
  gsl_spline_free(slfc_sp_ptr);
  gsl_interp_accel_free(slfc_acc_ptr);
 
}


// ---------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC STRUCTURE FACTOR
// ---------------------------------------------------------------------

void compute_dsf(double *SSn, double *phi_re, double *phi_im,
		 double GG, double *WW, input in){
    
  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double xx = in.dyn_xtarget;
  double xx2 = xx*xx;
  double ff1 = 4.0*lambda*in.rs/(M_PI*xx2);
  double ff2;
  double denom, denom_re, denom_im;
  
  for (int ii=0; ii<in.dyn_nW; ii++){

    if (xx == 0.0) {
      SSn[ii] = 0.0;
      continue;
    }
    
    if (WW[ii] == 0.0) {

      ff2 = 1.0/(1.0 + exp(xx2/(4.0*in.Theta) - in.mu));
      denom_re = 1.0 + ff1 * (1 - GG) * phi_re[ii];
      denom = denom_re*denom_re;
      SSn[ii] = in.Theta/(4.0*xx)*(ff2/denom);
      
    }
    else {
      
      ff2 = 1.0/(1.0 - exp(-WW[ii]/in.Theta));
      denom_re = 1.0 + ff1 * (1 - GG) * phi_re[ii];
      denom_im = ff1 * (1 - GG) * phi_im[ii];
      denom = denom_re*denom_re + denom_im*denom_im;
      SSn[ii] = (ff2/M_PI)*phi_im[ii]/denom;

    }

  }

}


// ---------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE INTERMEDIATE SCATTERING FUNCTION
// ---------------------------------------------------------------------

// Intermediate scattering function
void compute_isf(double *FF, double *tt, double *SSn,
		 double *WW, input in) {

  double dt = 1.0/ISF_NTAU; 
  double err;
  size_t nevals;

  // Declare accelerator and spline objects
  gsl_spline *dsf_sp_ptr;
  gsl_interp_accel *dsf_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  dsf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.dyn_nW);
  dsf_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(dsf_sp_ptr, WW, SSn, in.dyn_nW);

  // Integration workspace 
  gsl_integration_cquad_workspace *wsp 
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &isf;

  // Imaginary time
  for (int ii=0; ii<ISF_NTAU; ii++){
    tt[ii] = ii*dt;
  }
  
  // Compute intermediate scattering function
  for (int ii=0; ii<ISF_NTAU; ii++) {

    struct isf_params isfp = {in.Theta, tt[ii],
			      dsf_sp_ptr, dsf_acc_ptr};
    ff_int.params = &isfp;
    gsl_integration_cquad(&ff_int, 
			  WW[0], WW[in.dyn_nW-1], 
			  0.0, QUAD_REL_ERR, 
			  wsp, 
			  &FF[ii], &err, &nevals);

  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(dsf_sp_ptr);
  gsl_interp_accel_free(dsf_acc_ptr);
  
}


double isf(double WW, void *pp) {

  struct isf_params *params = (struct isf_params*)pp;
  double Theta = (params->Theta);
  double tau = (params->tau);
  gsl_spline *dsf_sp_ptr = (params->dsf_sp_ptr);
  gsl_interp_accel *dsf_acc_ptr = (params->dsf_acc_ptr);;
  double WW_T = WW/Theta;
  double WW_Tt = WW_T*tau;
  
  return 1.5*gsl_spline_eval(dsf_sp_ptr, WW, dsf_acc_ptr)*
    (exp(-WW_Tt) + exp(-WW_T + WW_Tt));
  
}

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

// write text files for output
void write_text_dynamic_stls(double *SSn, double *phi_re,
			     double *phi_im, double *WW,
			     input in){

  // Dynamic structure factor
  write_text_dsf(SSn, WW, in);

  // Intermediate scattering function
  write_text_isf(SSn, WW, in);

  // Ideal density response
  write_text_idr(phi_re, phi_im, WW, in);
  
}


// write static structure factor to text file
void write_text_dsf(double *SSn, double *WW, input in){

  FILE* fid;
  
  char out_name[100];
  sprintf(out_name, "dsf_rs%.3f_theta%.3f_x%.3f_%s.dat", in.rs, in.Theta,
	  in.dyn_xtarget, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the dynamic structure factor\n");
    exit(EXIT_FAILURE);
  }
    
  for (int ii = 0; ii<in.dyn_nW; ii++)
    fprintf(fid, "%.8e %.8e\n", WW[ii], SSn[ii]);
  
  fclose(fid);
  
}

// write intermediate scattering function to file
void write_text_isf(double *SSn, double *WW, input in){

  FILE* fid;
  double FF[ISF_NTAU];
  double tt[ISF_NTAU];
  
  char out_name[100];
  sprintf(out_name, "isf_rs%.3f_theta%.3f_x%.3f_%s.dat", in.rs, in.Theta,
	  in.dyn_xtarget, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the intermediate"
	    " scattering function\n");
    exit(EXIT_FAILURE);
  }

  compute_isf(FF, tt, SSn, WW, in);
  
  for (int ii = 0; ii<ISF_NTAU; ii++)
    fprintf(fid, "%.8e %.8e\n", tt[ii], FF[ii]);
  
  fclose(fid);
  
}


// write ideal density response to file
void write_text_idr(double *phi_re, double *phi_im, double *WW, input in){

  FILE* fid;

  char out_name[100];
  sprintf(out_name, "idr_rs%.3f_theta%.3f_x%.3f_%s.dat", in.rs, in.Theta,
	  in.dyn_xtarget, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  for (int ii = 0; ii < in.dyn_nW; ii++)
    fprintf(fid, "%.8e %.8e %.8e\n", WW[ii], phi_re[ii], phi_im[ii]);
  
  fclose(fid);
  
}


// write ideal density response to binary file
void write_bin_idr_2D(double *phi_re, double *phi_im,
		      input in){
  
  // Name of output file
  char out_name[100];
  sprintf(out_name, "dynamic_restart_rs%.3f_theta%.3f_%s.bin", in.rs, in.Theta,
	  in.theory);

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(out_name, "wb");
  if (fid == NULL) {
    fprintf(stderr,"Error while creating file to store the density"
	    " responses\n");
    exit(EXIT_FAILURE);
  }

  // Input data
  fwrite(&in.nx, sizeof(int), 1, fid);
  fwrite(&in.dx, sizeof(double), 1, fid);
  fwrite(&in.xmax, sizeof(double), 1, fid);
  fwrite(&in.dyn_nW, sizeof(int), 1, fid);
  fwrite(&in.dyn_dW, sizeof(double), 1, fid);
  fwrite(&in.dyn_Wmax, sizeof(double), 1, fid);
  fwrite(&in.Theta, sizeof(double), 1, fid);
  fwrite(&in.rs, sizeof(double), 1, fid);
  
  // Density response
  fwrite(phi_re, sizeof(double), in.nx * in.dyn_nW, fid);
  fwrite(phi_im, sizeof(double), in.nx * in.dyn_nW, fid);

  // Close binary file
  fclose(fid);

}

// read ideal density response from binary file
void read_bin_idr_2D(double *phi_re, double *phi_im, input in){
  
  // Variables
  size_t it_read;
  int nx_file;
  double dx_file;
  double xmax_file;
  int dyn_nW_file;
  double dW_file;
  double Wmax_file;
  double Theta_file;
  double rs_file;

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(in.dyn_restart_file, "rb");
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
  it_read += fread(&dyn_nW_file, sizeof(int), 1, fid);
  it_read += fread(&dW_file, sizeof(double), 1, fid);
  it_read += fread(&Wmax_file, sizeof(double), 1, fid);
  it_read += fread(&Theta_file, sizeof(double), 1, fid);
  it_read += fread(&rs_file, sizeof(double), 1, fid);
  check_bin_dynamic(dx_file, nx_file, xmax_file,
		    dW_file, dyn_nW_file, Wmax_file,
		    Theta_file, rs_file,
		    in, it_read, 8, fid, true, true, false);
  

  // Fixed component of the auxiliary density response 
  it_read += fread(phi_re, sizeof(double), nx_file * dyn_nW_file, fid);
  it_read += fread(phi_im, sizeof(double), nx_file * dyn_nW_file, fid);
  check_bin_dynamic(dx_file, nx_file, xmax_file,
		    dW_file, dyn_nW_file, Wmax_file,
		    Theta_file, rs_file,
		    in, it_read, 2*nx_file*dyn_nW_file + 8,
		    fid, false, true, false);
  
  
  // Close binary file
  fclose(fid);
	    
}


// Check consistency of the guess data
void check_bin_dynamic(double dx, int nx, double xmax,
		       double dW, int dyn_nW, double Wmax,
		       double Theta, double rs,
		       input in, size_t it_read,
		       size_t it_expected, FILE *fid,
		       bool check_grid, bool check_items,
		       bool check_eof){
  
  int buffer;
  
  // Check that the grid in the imported data is consistent with input
  if (check_grid) {

    if (nx != in.nx || fabs(dx-in.dx) > DBL_TOL || fabs(xmax-in.xmax) > DBL_TOL){
      fprintf(stderr,"Wave-vector grid from imported file is incompatible with input\n");
      fprintf(stderr,"Grid points (nx) : %d (input), %d (file)\n", in.nx, nx);
      fprintf(stderr,"Resolution (dx)  : %.16f (input), %.16f (file)\n", in.dx, dx);
      fprintf(stderr,"Cutoff (xmax)    : %.16f (input), %.16f (file)\n", in.xmax, xmax);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
    if (dyn_nW != in.dyn_nW || fabs(dW-in.dyn_dW) > DBL_TOL || fabs(Wmax-in.dyn_Wmax) > DBL_TOL){
      fprintf(stderr,"Frequency grid from imported file is incompatible with input\n");
      fprintf(stderr,"Grid points (dyn_nW) : %d (input), %d (file)\n", in.dyn_nW, dyn_nW);
      fprintf(stderr,"Resolution (dW)  : %.16f (input), %.16f (file)\n", in.dyn_dW, dW);
      fprintf(stderr,"Cutoff (Wmax)    : %.16f (input), %.16f (file)\n", in.dyn_Wmax, Wmax);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
    if (fabs(Theta-in.Theta) > DBL_TOL || fabs(rs-in.rs) > DBL_TOL){
      fprintf(stderr,"State point from imported file is incompatible with input\n");
      fprintf(stderr,"Degeneracy parameter (theta) : %f (input), %f (file)\n", in.Theta, Theta);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
    
  }

  // Check that all the expected items where read
  if (check_items) {
    if (it_read != it_expected ) {
      fprintf(stderr,"Error while reading file for the density response.\n");
      fprintf(stderr,"%ld Elements expected, %ld elements read\n", it_expected, it_read);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
  }
  
  // Check for end of file
  if (check_eof){
    it_read = fread(&buffer, sizeof(int), 1, fid); // Trigger end-of-file activation
    if (!feof(fid)) {
      fprintf(stderr,"Error while reading file for the density response.\n");
      fprintf(stderr,"Expected end of file, but there is still data left to read.\n");
      fclose(fid);
      exit(EXIT_FAILURE);
    }
  }
  
}

