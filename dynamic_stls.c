#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "solvers.h"
#include "utils.h"
#include "chemical_potential.h"
#include "stls.h"
#include "qstls.h"
#include "dynamic_stls.h"

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

  // Ideal density response
  if (verbose) printf("Normalized ideal Lindhard density calculation: ");
  compute_dynamic_idr(phi_re, phi_im, WW, in);
  if (verbose) printf("Done.\n");
  
  // Static local field correction
  if (verbose) printf("Static local field correction (from file): ");
  get_slfc(&GG, in);
  if (verbose) printf("Done.\n");
  
  // Dynamic structure factor
  if (verbose) printf("Dynamic structure factor calculation: ");
  compute_dsf(SSn, phi_re, phi_im, GG, WW, in);
  if (verbose) printf("Done.\n");
  
  // Output to file
  if (verbose) printf("Writing output files: ");
  write_text_dynamic_stls(SSn, WW, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_dynamic_stls_arrays(WW, phi_re, phi_im, SSn);

 
}

// -------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE SIZE OF THE FREQUENCY GRID
// -------------------------------------------------------------------

void get_frequency_grid_size(input *in){

  // Number of grid points in the frequency grid
  in->nW = (int)floor(in->dyn_Wmax/in->dyn_dW);
  
}
// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_stls_arrays(input in, double **WW, double **phi_re, 
			       double **phi_im, double **SSn){

  *WW = malloc( sizeof(double) * in.nW);
  if (*WW == NULL) {
    fprintf(stderr, "Failed to allocate memory for the frequency grid\n");
    exit(EXIT_FAILURE);
  }

  *phi_re = malloc( sizeof(double) * in.nW);
  if (*phi_re == NULL) {
    fprintf(stderr, "Failed to allocate memory for the real part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *phi_im = malloc( sizeof(double) * in.nW);
  if (*phi_im == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *SSn = malloc( sizeof(double) * in.nW);
  if (*SSn == NULL) {
    fprintf(stderr, "Failed to allocate memory for the dynamic structure factor\n");
    exit(EXIT_FAILURE);
  }

}

void free_dynamic_stls_arrays(double *WW, double *phi_re, 
			      double *phi_im, double *SSn){

  free(WW);
  free(phi_re);
  free(phi_im);
  free(SSn);
 
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
  for (int ii=1; ii < in->nW; ii++) WW[ii] = WW[ii-1] + in->dyn_dW;

}

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE IDEAL DENSITY RESPONSE
// ------------------------------------------------------------------

struct idr_params {

  double xx;
  double mu;
  double Theta;
  double WW;

};

// Ideal density response (real and imaginary part)
void compute_dynamic_idr(double *phi_re, double *phi_im,  double *WW,
			 input in) {

  // Real component
  compute_dynamic_idr_re(phi_re, WW, in);

  // Imaginary component
  compute_dynamic_idr_im(phi_im, WW, in);
  
}

// Real part of the ideal density response
void compute_dynamic_idr_re(double *phi_re, double *WW,
			     input in) {

  double xx = in.dyn_xtarget;
  double err;
  size_t nevals;
  
  // Integration workspace 
  gsl_integration_cquad_workspace *wsp 
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  if (WW == 0) ff_int.function = &idr_re_partial_x0;
  else ff_int.function = &idr_re_partial_xw;

  // Normalized ideal Lindhard density
  for (int ii=0; ii<in.nW; ii++) {

    struct idr_params phixwp = {xx, in.mu, in.Theta, WW[ii]};
    ff_int.params = &phixwp;
    gsl_integration_cquad(&ff_int, 
			  0.0, in.xmax, 
			  0.0, QUAD_REL_ERR, 
			  wsp, 
			  &phi_re[ii], &err, &nevals);

  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

// Imaginary part of the ideal density response
void compute_dynamic_idr_im(double *phi_im, double *WW,
			     input in) {

  double xx = in.dyn_xtarget;
  double ymin;
  double ymax;
  double err;
  size_t nevals;


  // Integration workspace 
  gsl_integration_cquad_workspace *wsp 
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  if (WW == 0) ff_int.function = &idr_im_partial_x0;
  else ff_int.function = &idr_im_partial_xw;

  // Normalized ideal Lindhard density
  for (int ii=0; ii<in.nW; ii++) {

    ymin = (xx/2.0) - WW[ii]/(2.0*xx);
    if (ymin < 0.0) ymin = -ymin;
    ymax = (xx/2.0) + WW[ii]/(2.0*xx);
    
    struct idr_params phixwp = {xx, in.mu, in.Theta, WW[ii]};
    ff_int.params = &phixwp;
    gsl_integration_cquad(&ff_int, 
			  ymin, ymax, 
			  0.0, QUAD_REL_ERR, 
			  wsp, 
			  &phi_im[ii], &err, &nevals);

  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

// Partial real part of the ideal density response (frequency = w, vector = x)
double idr_re_partial_xw(double yy, void *pp) {

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
double idr_re_partial_x0(double yy, void *pp) {

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
double idr_im_partial_xw(double yy, void *pp) {

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
double idr_im_partial_x0(double yy, void *pp) {
  return 0;
}


// ---------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE STATIC LOCAL FIELD CORRECTION (FROM FILE)
// ---------------------------------------------------------------------

void get_slfc(double *GG, input in){

  // Variables
  size_t it_read;
  char slfc_file_name[100];
  int nx_file;
  double dx_file;
  double xmax_file;
  double *SS_file = NULL;
  double *GG_file = NULL;
  double *xx_file = NULL;
  input in_file = in;
  gsl_spline *slfc_sp_ptr;
  gsl_interp_accel *slfc_acc_ptr;

  // Open binary file
  FILE *fid = NULL;
  if (strcmp(in.stls_guess_file, NO_FILE_STR)==0) {
    sprintf(slfc_file_name, "restart_rs%.3f_theta%.3f_%s.bin", in.rs, in.Theta, in.theory);
    fid = fopen(slfc_file_name, "rb");
  }
  else {
    fid = fopen(in.stls_guess_file, "rb");
  }
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file for the static local field correction\n");
    exit(EXIT_FAILURE);
  }

  // Initialize number of items read from input file
  it_read = 0;

  // Check that the data for the guess file is consistent
  it_read += fread(&nx_file, sizeof(int), 1, fid);
  it_read += fread(&dx_file, sizeof(double), 1, fid);
  it_read += fread(&xmax_file, sizeof(double), 1, fid);
  check_guess_stls(nx_file, dx_file, xmax_file, in, it_read, 3,
  		   fid, false, true, false);
  
  // Allocate temporary arrays to store the structural properties
  SS_file = malloc( sizeof(double) * nx_file);
  GG_file = malloc( sizeof(double) * nx_file);
  xx_file = malloc( sizeof(double) * nx_file);
  if (SS_file == NULL ||
      GG_file == NULL ||
      xx_file == NULL) {
    fprintf(stderr, "Failed to allocate memory for the data read"
  	    " from file\n");
    exit(EXIT_FAILURE);
  }
  
  // Static structure factor
  it_read += fread(SS_file, sizeof(double), nx_file, fid);

  // Static local field correction
  it_read += fread(GG_file, sizeof(double), nx_file, fid);

  // Check that all items where read and the end-of-file was reached
  check_guess_stls(nx_file, dx_file, xmax_file, in, it_read,
  		   2*nx_file + 3, fid, false, true, true);
 
  // Close binary file
  fclose(fid);

  // Wave-vector grid consistent with the data read from file
  in_file.nx = nx_file;
  in_file.dx = dx_file;
  in_file.xmax = xmax_file;
  wave_vector_grid(xx_file, &in_file);
    
  // Static local field correction for the wave-vector given in input
  slfc_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, nx_file);
  slfc_acc_ptr = gsl_interp_accel_alloc();
  gsl_spline_init(slfc_sp_ptr, xx_file, GG_file, nx_file);
  *GG = gsl_spline_eval(slfc_sp_ptr, in.dyn_xtarget, slfc_acc_ptr);
  
  // Free memory
  free(SS_file);
  free(GG_file);
  free(xx_file);
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
  
  for (int ii=0; ii<in.nW; ii++){

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

struct isf_params {

  double Theta;
  double tau;
  gsl_spline *dsf_sp_ptr;
  gsl_interp_accel *dsf_acc_ptr;
  
};


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
  dsf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nW);
  dsf_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(dsf_sp_ptr, WW, SSn, in.nW);

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
			  WW[0], WW[in.nW-1], 
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
void write_text_dynamic_stls(double *SSn, double *WW, input in){

  // Static structure factor
  write_text_dsf(SSn, WW, in);

  // Intermediate scattering function
  write_text_isf(SSn, WW, in);
  
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
  for (int ii = 0; ii < in.nW; ii++)
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
