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
#include "dynamic_qstls.h"

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

  // Ideal density response
  if (verbose) printf("Normalized ideal Lindhard density calculation: ");
  compute_dynamic_idr(phi_re, phi_im, WW, in);
  if (verbose) printf("Done.\n");
  
  // Static structure factor
  if (verbose) printf("Static structure factor (from file): ");
  get_ssf(&SS, &xx, &in);
  if (verbose) printf("Done.\n");

  // Auxiliary density response
  if (verbose) printf("Auxiliary density calculation: ");
  compute_dynamic_adr(psi_re, psi_im, WW, SS, xx, in);
  if (verbose) printf("Done.\n");

  // Dynamic structure factor
  if (verbose) printf("Dynamic structure factor calculation: ");
  compute_dsf_qstls(SSn, phi_re, phi_im, psi_re, psi_im, WW, in);
  if (verbose) printf("Done.\n");
  
  // Output to file
  if (verbose) printf("Writing output files: ");
  write_text_dynamic_qstls(SSn, WW, psi_re, psi_im, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_dynamic_stls_arrays(WW, phi_re, phi_im, SSn);
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
	    " the ideal density response\n");
    exit(EXIT_FAILURE);
  }
  
  *psi_im = malloc( sizeof(double) * in.nW);
  if (*psi_im == NULL) {
    fprintf(stderr, "Failed to allocate memory for the imaginary part of"
	    " the ideal density response\n");
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



// ---------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE STATIC STRUCTURE FACTOR (FROM FILE)
// ---------------------------------------------------------------------

void get_ssf(double **SS, double **xx, input *in){

  // Variables
  size_t it_read;
  char slfc_file_name[100];
  int nx_file;
  int nl_file;
  double dx_file;
  double xmax_file;
  input in_file = *in;

  // Open binary file
  FILE *fid = NULL;
  if (strcmp(in->stls_guess_file, NO_FILE_STR)==0) {
    sprintf(slfc_file_name, "restart_rs%.3f_theta%.3f_%s.bin", in->rs,
	    in->Theta, in->theory);
    fid = fopen(slfc_file_name, "rb");
  }
  else {
    fid = fopen(in->stls_guess_file, "rb");
  }
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file for the static local field correction\n");
    exit(EXIT_FAILURE);
  }

  // Initialize number of items read from input file
  it_read = 0;

  // Check that the data for the guess file is consistent
  it_read += fread(&nx_file, sizeof(int), 1, fid);
  it_read += fread(&nl_file, sizeof(int), 1, fid);
  it_read += fread(&dx_file, sizeof(double), 1, fid);
  it_read += fread(&xmax_file, sizeof(double), 1, fid);
  check_guess_qstls(nx_file, dx_file, xmax_file, nl_file, in->Theta,
		    *in, it_read, 4, fid, false, true, false);
  
  // Allocate temporary arrays to store the structural properties
  *SS = malloc( sizeof(double) * nx_file);
  *xx = malloc( sizeof(double) * nx_file);
  if (*SS == NULL ||
      *xx == NULL) {
    fprintf(stderr, "Failed to allocate memory for the data read"
  	    " from file\n");
    exit(EXIT_FAILURE);
  }
  
  // Static structure factor
  it_read += fread(*SS, sizeof(double), nx_file, fid);

  // Check that all items where read
  check_guess_qstls(nx_file, dx_file, xmax_file, nl_file, in->Theta,
		    *in, it_read, nx_file + 4, fid, false, true, false);
  
  // Close binary file
  fclose(fid);

  // Wave-vector grid consistent with the data read from file
  in_file.nx = nx_file;
  in_file.dx = dx_file;
  in_file.xmax = xmax_file;
  wave_vector_grid(*xx, &in_file);
  in->nx=in_file.nx;
  in->dx=in_file.dx;

}

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE AUXILIARY DENSITY RESPONSE
// ------------------------------------------------------------------

// Ideal density response (real and imaginary part)
void compute_dynamic_adr(double *psi_re, double *psi_im,
			 double *WW, double *SS,
			 double *xx, input in) {

  // Real component
  compute_dynamic_adr_re_part1(psi_re, WW, SS, xx, in);
  
  // Imaginary component
  compute_dynamic_adr_im_part1(psi_im, WW, SS, xx, in);
  
}


// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE REAL PART OF THE AUXILIARY DENSITY
// RESPONSE
// ------------------------------------------------------------------

// Real part of the auxiliary density response (part 1)

struct adr_re_part1_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *psi_re_part1_sp_ptr;
  gsl_interp_accel *psi_re_part1_acc_ptr;

};

void compute_dynamic_adr_re_part1(double *psi_re, double *WW,
				  double *SS, double *xx,
				  input in) {

  // Parallel calculations
  #pragma omp parallel
  {
  
    double err;
    size_t nevals;
    double *psi_re_part1  = malloc( sizeof(double) * in.nx);
    if (psi_re_part1 == NULL){
      fprintf(stderr, "Failed to allocate memory for calculation"
	      " of the real part of the auxiliary density"
	      " response function\n");
      exit(EXIT_FAILURE);
    }
    
    // Declare accelerator and spline objects
    gsl_spline *ssf_sp_ptr;
    gsl_interp_accel *ssf_acc_ptr;
    gsl_spline *psi_re_part1_sp_ptr;
    gsl_interp_accel *psi_re_part1_acc_ptr;
    
    // Allocate the accelerator and the spline objects
    ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    ssf_acc_ptr = gsl_interp_accel_alloc();
    psi_re_part1_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    psi_re_part1_acc_ptr = gsl_interp_accel_alloc();
    
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);
    
    // Loop over the frequency
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nW; ii++){
      
      // Integration function
      gsl_function ff_int_part1;
      ff_int_part1.function = &adr_re_part1_partial_xW;
      
      // Inner integrals
      compute_dynamic_adr_re_part2(psi_re_part1, WW[ii], xx, in);
    
      // Construct integrand
      gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);
      gsl_spline_init(psi_re_part1_sp_ptr, xx, psi_re_part1, in.nx);
      
      // Integral over w
      struct adr_re_part1_params ppart1 = {ssf_sp_ptr,
					   ssf_acc_ptr,
					   psi_re_part1_sp_ptr,
					   psi_re_part1_acc_ptr};
      ff_int_part1.params = &ppart1;
      gsl_integration_cquad(&ff_int_part1,
			    xx[0], xx[in.nx-1],
			    0.0, QUAD_REL_ERR,
			    wsp,
			    &psi_re[ii],
			    &err, &nevals);
      
    }
    
    // Free memory
    free(psi_re_part1);
    gsl_integration_cquad_workspace_free(wsp);
    gsl_spline_free(psi_re_part1_sp_ptr);
    gsl_interp_accel_free(psi_re_part1_acc_ptr);
    
  }
  
}

// Integrand for part 1 of the real auxiliary density response (vector = x, frequency = W)
double adr_re_part1_partial_xW(double ww, void* pp) {
  
  struct adr_re_part1_params* params = (struct adr_re_part1_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  gsl_spline* psi_re_part1_sp_ptr = (params->psi_re_part1_sp_ptr);
  gsl_interp_accel* psi_re_part1_acc_ptr = (params->psi_re_part1_acc_ptr);
  double ffp1 = gsl_spline_eval(psi_re_part1_sp_ptr, ww, psi_re_part1_acc_ptr);
  double ssfm1 = gsl_spline_eval(ssf_sp_ptr, ww, ssf_acc_ptr) - 1.0;

  return ww*ssfm1*ffp1;

}


// Real part of the auxiliary density response (part 2)

struct adr_re_part2_params {

  double mu;
  double Theta;
  gsl_spline *psi_re_part2_sp_ptr;
  gsl_interp_accel *psi_re_part2_acc_ptr;

};

void compute_dynamic_adr_re_part2(double *psi_re_part1, double WW,
				  double *ww, input in) {

  double err;
  size_t nevals;
  double *psi_re_part2  = malloc( sizeof(double) * in.nx);
  if (psi_re_part2 == NULL){
    fprintf(stderr, "Failed to allocate memory for calculation"
	    " of the real part of the auxiliary density"
	    " response function\n");
    exit(EXIT_FAILURE);
  }
  
  // Declare accelerator and spline objects
  gsl_spline *psi_re_part2_sp_ptr;
  gsl_interp_accel *psi_re_part2_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  psi_re_part2_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  psi_re_part2_acc_ptr = gsl_interp_accel_alloc();
    
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
    
 
  // Integration function
  gsl_function ff_int_part2;
  if (WW == 0){
    ff_int_part2.function = &adr_re_part2_partial_xw0;
  }
  else {
    ff_int_part2.function = &adr_re_part2_partial_xwW;
  }
  
  // Loop over w (wave-vector)
  for (int ii=0; ii<in.nx; ii++){

    // Inner integral
    compute_dynamic_adr_re_part3(psi_re_part2, WW, ww[ii], ww, in);
	
    // Construct integrand
    gsl_spline_init(psi_re_part2_sp_ptr, ww, psi_re_part2, in.nx);
    
    // Integrate
    struct adr_re_part2_params ppart2 = {in.mu,in.Theta,
					 psi_re_part2_sp_ptr,
					 psi_re_part2_acc_ptr};
    ff_int_part2.params = &ppart2;
    gsl_integration_cquad(&ff_int_part2,
			  ww[0], ww[in.nx-1],
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &psi_re_part1[ii],
			  &err, &nevals);
  }

  // Free memory
  free(psi_re_part2);
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(psi_re_part2_sp_ptr);
  gsl_interp_accel_free(psi_re_part2_acc_ptr);
  
}

// Integrand for part 2 of the real auxiliary density response (vectors = {x,w}, frequency = W)
double adr_re_part2_partial_xwW(double qq, void* pp) {
  
  struct adr_re_part2_params* params = (struct adr_re_part2_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  gsl_spline* psi_re_part2_sp_ptr = (params->psi_re_part2_sp_ptr);
  gsl_interp_accel* psi_re_part2_acc_ptr = (params->psi_re_part2_acc_ptr);
  double qq2 = qq*qq;
  double ffp2 = gsl_spline_eval(psi_re_part2_sp_ptr, qq, psi_re_part2_acc_ptr);

  return -(3.0/8.0)*qq/(exp(qq2/Theta - mu) + 1.0)*ffp2;

}

// Integrand for part 2 of the real auxiliary density response (vectors = {x,w}, frequency = 0)
double adr_re_part2_partial_xw0(double qq, void* pp) {

  struct adr_re_part2_params* params = (struct adr_re_part2_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  gsl_spline* psi_re_part2_sp_ptr = (params->psi_re_part2_sp_ptr);
  gsl_interp_accel* psi_re_part2_acc_ptr = (params->psi_re_part2_acc_ptr);
  double qq2 = qq*qq;
  double ffp2 = gsl_spline_eval(psi_re_part2_sp_ptr, qq, psi_re_part2_acc_ptr);

  return  -(3.0/(4.0*Theta))*qq/(exp(qq2/Theta - mu)
				  + exp(-qq2/Theta + mu) + 2.0)*ffp2;

}

 
// Real part of the auxiliary density response (part 3)

struct adr_re_part3_params {

  double ww;
  double xx;
  double WW;
  double qq;

};

void compute_dynamic_adr_re_part3(double *psi_re_part2, double WW,
				  double ww, double *qq, input in) {

  double err;
  size_t nevals;

  // Integration limits
  double xx = in.dyn_xtarget;
  double xx2 = xx*xx;
  double xw = xx*ww;
  double tmin = xx2 - xw;
  double tmax = xx2 + xw;
    
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
     
  // Integration function
  gsl_function ff_int_part3;
  if (WW == 0){
    ff_int_part3.function = &adr_re_part3_partial_xwq0;
  }
  else {
    ff_int_part3.function = &adr_re_part3_partial_xwqW;
  }
  
    
  // Loop over q (wave-vector)
  for (int ii=0; ii<in.nx; ii++) {
    
    if (qq[ii] > 0.0){
      
      // Integration over t
      struct adr_re_part3_params ppart3 = {ww,xx,WW,qq[ii]};
      ff_int_part3.params = &ppart3;
      gsl_integration_cquad(&ff_int_part3,
			    tmin, tmax,
			    0.0, QUAD_REL_ERR,
			    wsp,
			    &psi_re_part2[ii],
			    &err, &nevals);
    }
    
    else psi_re_part2[ii] = 0.0;
    
  }

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}


// Integrand for part 2 of the real auxiliary density response (vectors = {x,q,w}, frequency = W)
double adr_re_part3_partial_xwqW(double tt, void* pp) {
  
  struct adr_re_part3_params* params = (struct adr_re_part3_params*)pp;
  double xx = (params->xx);
  double qq = (params->qq);
  double ww = (params->ww);
  double WW = (params->WW);
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double WW2 = WW*WW; 
  double txq = 2.0*xx*qq; 
  double txqpt = txq + tt;
  double txqmt = txq - tt;
  double txqpt2 = txqpt*txqpt;
  double txqmt2 = txqmt*txqmt;
  double logarg = (txqpt2 - WW2)/(txqmt2 - WW2);

  if (logarg < 0.0) logarg = -logarg;
  if (xx == 0 || qq == 0)
    return 0;
  else
    return 1.0/(2.0*tt + ww2 - xx2)*log(logarg);

}

// Integrand for part 2 of the real auxiliary density response (vectors = {x,q,w}, frequency = 0)
double adr_re_part3_partial_xwq0(double tt, void* pp) {

  struct adr_re_part3_params* params = (struct adr_re_part3_params*)pp;
  double xx = (params->xx);
  double qq = (params->qq);
  double ww = (params->ww);
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double qq2 = qq*qq; 
  double tt2 = tt*tt;
  double txq = 2.0*xx*qq;
  double logarg;

  if (xx == 0 || qq == 0){
    return 0;
  }
  else if  (tt == txq){
    return 2.0*qq2/(ww2 + 2.0*txq - xx2);
  }
  else {
    
    logarg = (tt + txq)/(tt - txq);
    if (logarg < 0.0) logarg = -logarg;
    return 1.0/(2.0*tt + ww2 - xx2)*((qq2 - tt2/(4.0*xx2))
    				     *log(logarg) + qq*tt/xx);
				     
  }

}


// ------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE IMAGINARY PART OF THE AUXILIARY 
// DENSITY RESPONSE
// ------------------------------------------------------------------

// Imaginary part of the auxiliary density response (part 1)

struct adr_im_part1_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *psi_im_part1_sp_ptr;
  gsl_interp_accel *psi_im_part1_acc_ptr;

};

void compute_dynamic_adr_im_part1(double *psi_im, double *WW,
				  double *SS, double *xx,
				  input in) {

  // Parallel calculations
  #pragma omp parallel
  {
  
    double err;
    size_t nevals;
    double *psi_im_part1  = malloc( sizeof(double) * in.nx);
    if (psi_im_part1 == NULL){
      fprintf(stderr, "Failed to allocate memory for calculation"
	      " of the imaginary part of the auxiliary density"
	      " response function\n");
      exit(EXIT_FAILURE);
    }
    
    // Declare accelerator and spline objects
    gsl_spline *ssf_sp_ptr;
    gsl_interp_accel *ssf_acc_ptr;
    gsl_spline *psi_im_part1_sp_ptr;
    gsl_interp_accel *psi_im_part1_acc_ptr;
    
    // Allocate the accelerator and the spline objects
    ssf_sp_ptr = gsl_spline_alloc(gsl_interp_linear, in.nx);
    ssf_acc_ptr = gsl_interp_accel_alloc();
    psi_im_part1_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    psi_im_part1_acc_ptr = gsl_interp_accel_alloc();
    
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);
    
    // Loop over the frequency
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nW; ii++){

      if (WW[ii] == 0.0) {
	psi_im[ii] = 0.0;
	continue;
      }
	  
      // Integration function
      gsl_function ff_int_part1;
      ff_int_part1.function = &adr_im_part1_partial_xW;
      
      // Inner integrals
      compute_dynamic_adr_im_part2(psi_im_part1, WW[ii], xx, in);
    
      // Construct integrand
      gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);
      gsl_spline_init(psi_im_part1_sp_ptr, xx, psi_im_part1, in.nx);
      
      // Integral over w
      struct adr_im_part1_params ppart1 = {ssf_sp_ptr,
					   ssf_acc_ptr,
					   psi_im_part1_sp_ptr,
					   psi_im_part1_acc_ptr};
      ff_int_part1.params = &ppart1;
      gsl_integration_cquad(&ff_int_part1,
			    xx[0], xx[in.nx-1],
			    0.0, QUAD_REL_ERR,
			    wsp,
			    &psi_im[ii],
			    &err, &nevals);
      
    }
    
    // Free memory
    free(psi_im_part1);
    gsl_integration_cquad_workspace_free(wsp);
    gsl_spline_free(psi_im_part1_sp_ptr);
    gsl_interp_accel_free(psi_im_part1_acc_ptr);
    
  }
  
}

// Integrand for part 1 of the imaginary auxiliary density response (vector = x, frequency = W)
double adr_im_part1_partial_xW(double ww, void* pp) {
  
  struct adr_im_part1_params* params = (struct adr_im_part1_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  gsl_spline* psi_im_part1_sp_ptr = (params->psi_im_part1_sp_ptr);
  gsl_interp_accel* psi_im_part1_acc_ptr = (params->psi_im_part1_acc_ptr);
  double ffp1 = gsl_spline_eval(psi_im_part1_sp_ptr, ww, psi_im_part1_acc_ptr);
  double ssfm1 = gsl_spline_eval(ssf_sp_ptr, ww, ssf_acc_ptr) - 1.0;

  return ww*ssfm1*ffp1;

}


// Imaginary part of the auxiliary density response (part 2)

struct adr_im_part2_params {

  double ww;
  double xx;
  gsl_spline *psi_im_part2_sp_ptr;
  gsl_interp_accel *psi_im_part2_acc_ptr;
  
};

void compute_dynamic_adr_im_part2(double *psi_im_part1, double WW,
				  double *ww, input in) {

  double err;
  size_t nevals;

  // Integration limits
  double xx = in.dyn_xtarget;
  double uu[ADR_IM_NU];
  double du = 2.0/(ADR_IM_NU - 1);
  double psi_im_part2[ADR_IM_NU];
  
  // Declare accelerator and spline objects
  gsl_spline *psi_im_part2_sp_ptr;
  gsl_interp_accel *psi_im_part2_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  psi_im_part2_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, ADR_IM_NU);
  psi_im_part2_acc_ptr = gsl_interp_accel_alloc();
  
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
     
  // Integration function
  gsl_function ff_int_part2;
  ff_int_part2.function = &adr_im_part2_partial_xwW;
  

  // Fill array with integration variable (u)
  for (int ii=0; ii<ADR_IM_NU; ii++){
    uu[ii] = -1 + du*ii;
  }
  
  // Loop over w (wave-vector)
  for (int ii=0; ii<in.nx; ii++) {

    // Inner integral
    compute_dynamic_adr_im_part3(psi_im_part2, WW, ww[ii], ww, uu, in);

    // Construct integrand
    gsl_spline_init(psi_im_part2_sp_ptr, uu, psi_im_part2, ADR_IM_NU);
    
    // Integration over u (wave-vector squared)
    struct adr_im_part2_params ppart2 = {ww[ii], xx,
					 psi_im_part2_sp_ptr,
                                         psi_im_part2_acc_ptr};
    
    ff_int_part2.params = &ppart2;
    gsl_integration_cquad(&ff_int_part2,
			  uu[0], uu[ADR_IM_NU-1],
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &psi_im_part1[ii],
			  &err, &nevals);
  }
 
    

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(psi_im_part2_sp_ptr);
  gsl_interp_accel_free(psi_im_part2_acc_ptr);
  
}


// Integrand for part 2 of the imaginary auxiliary density response (vectors = {x,q,w}, frequency = W)
double adr_im_part2_partial_xwW(double uu, void* pp) {
  
  struct adr_im_part2_params* params = (struct adr_im_part2_params*)pp;
  double xx = (params->xx);
  double ww = (params->ww);
  gsl_spline* psi_im_part2_sp_ptr = (params->psi_im_part2_sp_ptr);
  gsl_interp_accel* psi_im_part2_acc_ptr = (params->psi_im_part2_acc_ptr);  
  double xx2 = xx*xx;
  double ww2 = ww*ww;
  double denom = xx2 +  ww2 - 2.0*xx*ww*uu;
  double ffp2 = gsl_spline_eval(psi_im_part2_sp_ptr, uu, psi_im_part2_acc_ptr);

  return (3.0*M_PI/8)*xx*ww*ffp2/denom;
  
}


// Imaginary part of the auxiliary density response (part 3)

struct adr_im_part3_params {

  double mu;
  double Theta;
  double xx;
  double ww;
  double uu;
  double WW;
  
};

void compute_dynamic_adr_im_part3(double *psi_im_part2, double WW,
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
  gsl_function ff_int_part3;
  ff_int_part3.function = &adr_im_part3_partial_xwuW;
  
  // Loop over u (wave-vector squared)
  for (int ii=0; ii<ADR_IM_NU; ii++){
    
    // Integration limits
    x2mxwu = xx2 - xw*uu[ii];
    if (x2mxwu < 0.0) x2mxwu = -x2mxwu;
    q_min = (WW - x2mxwu)/(2.0*xx);
    if (q_min < 0.0) q_min = -q_min;
    q_max = (WW + x2mxwu)/(2.0*xx);
    
    // Integrate over q (wave-vector)
    struct adr_im_part3_params ppart3 = {in.mu,in.Theta, xx, ww, uu[ii], WW};
    ff_int_part3.params = &ppart3;
    gsl_integration_cquad(&ff_int_part3,
			  q_min, q_max,
			  0.0, QUAD_REL_ERR,
			  wsp,
			  &psi_im_part2[ii],
			  &err, &nevals);
  }

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

// Integrand for part 3 of the imaginary auxiliary density response (vectors = {x,w}, frequency = W)
double adr_im_part3_partial_xwuW(double qq, void* pp) {
  
  struct adr_im_part3_params* params = (struct adr_im_part3_params*)pp;
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
  int out1 = 0;
  int out2 = 0;
  
  if (qq2 > hh1) 
    out1 = 1;

  if (qq2 > hh2)
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

    ff2 = 1.0/(1.0 - exp(-WW[ii]/in.Theta));
    numer = phi_im[ii] + ff1*(phi_re[ii]*psi_im[ii] -
    			      phi_im[ii]*psi_re[ii]);
    denom_re = 1.0 + ff1 * (phi_re[ii] - psi_re[ii]);
    denom_im = ff1 * (phi_im[ii] - psi_im[ii]);
    denom = denom_re*denom_re + denom_im*denom_im;

    if (xx == 0.0)
      SSn[ii] = 0.0;
    else
      SSn[ii] = (ff2/M_PI)*numer/denom;

  }

}


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

// write text files for output
void write_text_dynamic_qstls(double *SSn, double *WW, double *psi_re,
			      double *psi_im, input in){

  // Static structure factor
  write_text_dsf(SSn, WW, in);

  FILE* fid;
  
  char out_name[100];
  sprintf(out_name, "psire_rs%.3f_theta%.3f_x%.3f_%s.dat", in.rs, in.Theta,
	  in.dyn_xtarget, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the dynamic structure factor\n");
    exit(EXIT_FAILURE);
  }
  for (int ii = 0; ii < in.nW; ii++)
    fprintf(fid, "%.8e %.8e\n", WW[ii], psi_re[ii]);
  
  fclose(fid);

  
  sprintf(out_name, "psiim_rs%.3f_theta%.3f_x%.3f_%s.dat", in.rs, in.Theta,
	  in.dyn_xtarget, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the dynamic structure factor\n");
    exit(EXIT_FAILURE);
  }
  for (int ii = 0; ii < in.nW; ii++)
    fprintf(fid, "%.8e %.8e\n", WW[ii], psi_im[ii]);
  
  fclose(fid);
  
}

