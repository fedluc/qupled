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

  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
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
  if (verbose) printf("Internal energy: %f\n",compute_internal_energy(SS, in));
  
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

void compute_slfc_hnc(double *GG_new, double *GG, 
		      double *SS, double *xx, input in) {

  double wmax, wmin;
  double *GGu  = malloc( sizeof(double) * in.nx);

  // Static local field correction
  for (int ii=0; ii<in.nx; ii++) {
   
    

  }

}

double slfc_u(double uu, void* pp) {

    if (uu >= u_min_cut && uu <= u_max_cut && uu > 0.0)
      return (1.0/uu) * gsl_spline_eval(GGu_sp_ptr, uu, GGu_acc_ptr)
	*(1 - (gsl_spline_eval(ssf_sp_ptr, uu, ssf_acc_ptr) - 1.0)
	  *(gsl_spline_eval(slfc_sp_ptr, uu, slfc_acc_ptr) - 1.0));
    else 
      return 0;
}

double slfc_w(double ww, void* pp) {

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

