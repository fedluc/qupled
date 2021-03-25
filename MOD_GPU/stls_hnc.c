#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "stls.h"
#include "stls_hnc.h"


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS-HNC EQUATIONS
// -------------------------------------------------------------------

void solve_stls_hnc(input in, bool verbose) {

  // Arrays for STLS solution
  float *xx = NULL; 
  float *phi = NULL;
  float *GG = NULL;
  float *GG_new = NULL;
  float *SS = NULL;
  float *SSHF = NULL;

  // Solve STLS equation for initial guess
  if (verbose) printf("Solution of classical STLS for initial guess:\n");
  float a_mix_hold = in.a_mix;
  in.a_mix = 0.1;
  solve_stls(in, false, &xx, &SS, &SSHF, &GG, &GG_new, &phi);
  in.a_mix = a_mix_hold;
  if (verbose) printf("Done.\n");
  for (int ii=0; ii<in.nx; ii++){
    GG[ii] = 0.0;
  }

  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  float iter_err = 1.0;
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
      printf("Elapsed time: %f seconds\n", ((float)toc - (float)tic) / CLOCKS_PER_SEC);
      printf("Residual error: %.5e\n", iter_err);
      fflush(stdout);
    }
  }
  if (verbose) printf("Done.\n");
  
  // Internal energy
  if (verbose) printf("Internal energy: %f\n",compute_uex(SS, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text_hnc(SS, GG, xx, in);
  if (verbose) printf("Done.\n");

  // Output to variable or free memory
  free_stls_arrays(xx, true, phi, true, GG, true, 
		   GG_new, true, SS, true, SSHF, true);
 
 
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_hnc(float *GG_new, float *GG, 
		      float *SS, float *xx, input in) {

  float kmax, kmin, fu, 
    xx2, uu, uu2, ww2, ww;

  // Static local field correction
  for (int ii=0; ii<in.nx; ii++) {
   
    xx2 = xx[ii]*xx[ii];
    GG_new[ii] = 0.0;
    
    for (int jj=0; jj<in.nx; jj++){

      uu = xx[jj];
      uu2 = uu*uu;
      kmin = ii-jj;
      if (kmin < 0) kmin = -kmin;
      kmax = ii+jj;
      if (kmax > in.nx) kmax = in.nx;
      
      fu = 0.0;
      for (int kk=kmin; kk<kmax; kk++){

	ww = xx[kk];
	ww2 = ww*ww;
	fu += (ww2 - uu2 - xx2) * ww * (SS[kk] - 1.0);
	
      }
      
      GG_new[ii] += (1.0 - (GG[jj] - 1.0) * (SS[jj] - 1.0))*in.dx*fu/uu;
      
    }
    
    GG_new[ii] *= 3.0*in.dx/(8.0*xx[ii]);
    
  }

}


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------


// write text files with SSF and SLFC
void write_text_hnc(float *SS, float *GG, float *xx, input in){


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

