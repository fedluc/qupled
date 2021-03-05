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
  double a_mix_hold = in.a_mix;
  in.a_mix = 0.1;
  solve_stls(in, false, &xx, &SS, &SSHF, &GG, &GG_new, &phi);
  in.a_mix = a_mix_hold;
  if (verbose) printf("Done.\n");

  /* // Initial guess for Static structure factor (SSF) and static-local field correction (SLFC) */
  /* for (int ii=0; ii < in.nx; ii++) { */
  /*   GG[ii] = 0.0; */
  /*   GG_new[ii] = 0.0; */
  /* } */
  /* compute_ssf(SS, SSHF, GG, phi, xx, in); */

  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min ) {
    
    // Start timing
    clock_t tic = clock();
    
    // Update SLFC
    compute_slfc_hnc(GG_new, GG, SS, xx, in);
    
    // Update diagnostic
    iter_err = 0.0;
    iter_counter++;
    for (int ii=0; ii<in.nx; ii++) {
      iter_err += (GG_new[ii] - GG[ii]) * (GG_new[ii] - GG[ii]);
      GG[ii] = in.a_mix*GG_new[ii] + (1-in.a_mix)*GG[ii];
    }
    iter_err = sqrt(iter_err);
    
    // Update SSF
    compute_ssf(SS, SSHF, GG, phi, xx, in);
    
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
  write_text_hnc(SS, GG, xx, in);
  if (verbose) printf("Done.\n");

  // Output to variable or free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
 
 
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
  double u_min_cut;
  double u_max_cut;

};

struct slfcw_params {

  double xx;
  double uu;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  double w_min_cut;
  double w_max_cut;

};


void compute_slfc_hnc(double *GG_new, double *GG, 
		      double *SS, double *xx, input in) {

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
	
	if (xx[jj] >  0.0) {

	  struct slfcw_params slfcwp = {xx[ii], xx[jj], 
					ssf_sp_ptr, ssf_acc_ptr,
					xx[0], xx[in.nx-1]};
	  wmin = xx[jj] - xx[ii];
	  if (wmin < 0.0) wmin = -wmin;
	  wmax = GSL_MIN(xx[in.nx-1], xx[ii]+xx[jj]);
	  fw_int.params = &slfcwp;
	  gsl_integration_cquad(&fw_int,
				wmin, wmax,
				0.0, 1e-6,
				wsp,
				&GGu[jj], &err, &nevals);
	}

	else GGu[jj] = 0.0; 

      }

      // Interpolate result of integration over w
      gsl_spline_init(GGu_sp_ptr, xx, GGu, in.nx);    
      
      // Evaluate integral over u
      struct slfcu_params slfcup = {ssf_sp_ptr, ssf_acc_ptr,
				    slfc_sp_ptr, slfc_acc_ptr,
				    GGu_sp_ptr, GGu_acc_ptr,
				    xx[0], xx[in.nx-1]};
      fu_int.params = &slfcup;
      gsl_integration_cquad(&fu_int,
			    xx[0], xx[in.nx-1],
			    0.0, 1e-6,
			    wsp,
			    &GG_new[ii], &err, &nevals);
      
      GG_new[ii] *= 3.0/(8.0*xx[ii]);

    }

    else GG_new[ii] = 0.0;

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

double slfc_u(double uu, void* pp) {

    struct slfcu_params* params = (struct slfcu_params*)pp;
    gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
    gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
    gsl_spline* slfc_sp_ptr = (params->slfc_sp_ptr);
    gsl_interp_accel* slfc_acc_ptr = (params->slfc_acc_ptr);
    gsl_spline* GGu_sp_ptr = (params->GGu_sp_ptr);
    gsl_interp_accel* GGu_acc_ptr = (params->GGu_acc_ptr);
    double u_min_cut = (params->u_min_cut);
    double u_max_cut = (params->u_max_cut);

    if (uu >= u_min_cut && uu <= u_max_cut)
      return (1.0/uu) * gsl_spline_eval(GGu_sp_ptr, uu, GGu_acc_ptr)
	*(1 - (gsl_spline_eval(ssf_sp_ptr, uu, ssf_acc_ptr) - 1.0)
	  *(gsl_spline_eval(slfc_sp_ptr, uu, slfc_acc_ptr) - 1.0));
    else 
      return 0;
}

double slfc_w(double ww, void* pp) {

    struct slfcw_params* params = (struct slfcw_params*)pp;
    double xx = (params->xx);
    double uu = (params->uu);
    double w_min_cut = (params->w_min_cut);
    double w_max_cut = (params->w_max_cut);
    gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
    gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
    double ww2 = ww*ww, xx2 = xx*xx, uu2 = uu*uu;
    
    if (ww >= w_min_cut && ww <= w_max_cut)
      return (ww2 - uu2 - xx2)*ww
	*(gsl_spline_eval(ssf_sp_ptr, ww, ssf_acc_ptr) - 1.0);
    else 
      return 0;
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

