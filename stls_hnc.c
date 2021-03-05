#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "stls.h"
#include "stls_hnc.h"


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS-HNC EQUATIONS
// -------------------------------------------------------------------

void solve_stls_hnc(input in, bool verbose) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *GG_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;

  // Solve STLS equation for initial guess
  if (verbose) printf("Solution of classical STLS for initial guess:\n");
  solve_stls(in, false, &xx, &SS, &SSHF, &GG, &GG_new, &phi); 
  if (verbose) printf("Done.\n");
  
  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min ) {
    
    // Start timing
    clock_t tic = clock();
    
    // Update SLFC
    compute_slfc_hnc(GG_new, SS, xx, in);
    
    // Update diagnostic
    iter_err = 0.0;
    iter_counter++;
    for (int ii=0; ii<in.nx; ii++) {
      iter_err += (GG_new[ii] - GG[ii]) * (GG_new[ii] - GG[ii]);
      GG[ii] = in.a_mix*GG_new[ii] + (1-in.a_mix)*GG[ii];
    }
    iter_err = sqrt(iter_err);
    
    // Update SSF
    compute_ssf(SS, SSHF, AA, GG, phi, xx, in);
    
    // End timing
    clock_t toc = clock();
    
    // Print diagnostic
    if (verbose) {
      printf("--- iteration %d ---\n", iter_counter);
      printf("Elapsed time: %f seconds\n", ((double)toc - (double)tic) / CLOCKS_PER_SEC);
      printf("Residual error: %.5e\n", iter_err);
    }
  }
  if (verbose) printf("Done.\n");
  
  // Internal energy
  if (verbose) printf("Internal energy: %f\n",compute_internal_energy(SS, xx, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text(SS, GG, xx, in);
  if (init_flag) write_bin(phi, SSHF, AA, in);
  if (verbose) printf("Done.\n");

  // Output to variable or free memory
  if (xx_out != NULL) {
    xx_out = xx;
    SS_out = SS;
    SSHF_out = SSHF;
    GG_out = GG;
    phi_out = phi;
  }
  else{
    free_stls_arrays(xx, phi, AA, GG, GG_new, SS, SSHF);
  }
 
 
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

struct slfcu_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *slfc_sp_ptr;
  gsl_interp_accel *slfc_acc_ptr;
  gsl_spline *GGu_sp_ptr;
  gsl_interp_accel *GGu_acc_ptr;
  

};

struct slfcw_params {

  double xx;
  double uu;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};


void compute_slfc_hnc(double *GG_new, double *GG, double *SS, double *xx, input in) {

  double err;
  size_t nevals;
  double wmax, wmin;
  double *GGu  = malloc( sizeof(double) * in.nx);

  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *slfc_sp_ptr;
  gsl_interp_accel *slfc_acc_ptr;
  gsl_spline *GGu_sp_ptr;
  gsl_interp_accel *GGu_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  slfc_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  slfc_acc_ptr = gsl_interp_accel_alloc();
  GGu_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  GGu_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);    
  gsl_spline_init(slfc_sp_ptr, xx, GG, in.nx);    

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(in.nx);

  // Integration function
  gsl_function fu_int, fw_int;
  fu_int.function = &slfc_u;
  fw_int.function = &slfc_w;

  // Static local field correction
  // Integration over u
  for (int ii=0; ii<in.nx; ii++) {
    
    
    if (xx[ii] > 0.0){
      
      // Integration over w
      for (int jj=0; jj<in.nx; jj++){
 
	if (xx[jj] == 0) {

	  struct slfcw_params slfcwp = {xx[ii], xx[jj], ssf_sp_ptr, ssf_acc_ptr};
	  wmin = xx[jj] - xx[ii];
	  if (wmin < 0.0) wmin = -wmin;
	  wmax = GSL_MIN(xx[in.nx-1], xx[jj]+xx[ii]);
	  fw_int.params = &slfcwp;
	  gsl_integration_cquad(&fw_int,
				wmin, wmax,
				0.0, 1e-6,
				wsp,
				&GGu[jj], &err, &nevals);

	}
	else GGu[jj] = 0.0; 

	// Interpolate result of integration over w
	gsl_spline_init(GGu_sp_ptr, xx, GGu, in.nx);    
	  
	// Evaluate integral over u
	struct slfcu_params slfcup = {ssf_sp_ptr, ssf_acc_ptr,
				      sflc_sp_ptr, slfc_acc_ptr,
				      GGu_sp_ptr, GGu_sp_ptr};
	fu_int.params = &slfcup;
	gsl_integration_cquad(&fu_int,
			      xx[0], xx[in.nx - 1],
			      0.0, 1e-6,
			      wsp,
			      &GG_new[ii], &err, &nevals);
	
	GG[ii] *= -3.0/(8.0*xx[ii]);
      }

    }

    else GG[ii] = 0.0;

  }

  // Free memory
  free(GGu);
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(ssf_sp_ptr);
  gsl_interp_accel_free(ssf_acc_ptr);
  gsl_spline_free(slfc_sp_ptr);
  gsl_interp_accel_free(slfc_acc_ptr);
  gsl_spline_free(GGu_sp_ptr);
  gsl_interp_accel_free(GGu_acc_ptr);

}

double slfc_u(double yy, void* pp) {

    struct slfc_params* params = (struct slfc_params*)pp;
    double xx = (params->xx);
    gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
    gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
    double yy2 = yy * yy, xx2 = xx * xx;

    if (xx > yy){
      return yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0) 
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(xx - yy)));
    } 
    else if (xx < yy) {
      return yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0) 
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(yy - xx)));
    }
    else {
      	return 0;
    }

}

double slfc_w(double ww, void* pp) {

    struct slfc_params* params = (struct slfc_params*)pp;
    double xx = (params->xx);
    double uu = (params->uu);
    gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
    gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
    double ww2 = ww * ww, xx2 = xx * xx, uu2 = uu*uu;
    
    return (ww2 - uu2 - xx2)*ww
      *gsl_spline_eval(ssf_sp_ptr, ww, ssf_acc_ptr) - 1.0)

}


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------


// write text files with SSF and SLFC
void write_text_hnc(double *SS, double *GG, double *xx, input in){


    FILE* fid;
    
    // Output for SSF
    fid = fopen("ssf_STLS_HNC.dat", "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static structure factor");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nx; ii++)
    {
        fprintf(fid, "%.8e %.8e\n", xx[ii], SS[ii]);
    }
    fclose(fid);

    // Output for SLFC
    fid = fopen("slfc_STLS_HNC.dat", "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static local field correction");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nx; ii++)
    {
        fprintf(fid, "%.8e %.8e\n", xx[ii], GG[ii]);
    }
    fclose(fid);

}

