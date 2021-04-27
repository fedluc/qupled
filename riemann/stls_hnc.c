#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include "solvers.h"
#include "stls.h"
#include "stls_hnc.h"

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS-HNC EQUATIONS
// -------------------------------------------------------------------

void solve_stls_hnc(input in, bool verbose, bool iet) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *GG_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;
  double *bf = NULL;

  // Allocate arrays
  alloc_stls_arrays(in, &xx, &phi, &GG, &GG_new, &SS, &SSHF);

  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_stls_arrays(&in, xx, phi, SSHF, verbose);
  bf = malloc( sizeof(double) * in.nx);
  compute_bf(bf, xx, in, iet);

  // Initial guess for Static structure factor (SSF) and static-local field correction (SLFC)
  if (strcmp(in.guess_file,"NO_FILE")==0){
    for (int ii=0; ii < in.nx; ii++) {
      GG[ii] = 0.0;
      GG_new[ii] = 1.0;
    }
    compute_ssf(SS, SSHF, GG, phi, xx, in);
  }
  else {
    read_guess(SS, GG, in);
  }
   
  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    double tic = omp_get_wtime();
    
    // Update SSF
    compute_ssf(SS, SSHF, GG, phi, xx, in);
    
    // Update SLFC
    compute_slfc_hnc(GG_new, GG, SS, bf, xx, in);
    
    // Update diagnostic
    iter_err = 0.0;
    iter_counter++;
    for (int ii=0; ii<in.nx; ii++) {
      iter_err += (GG_new[ii] - GG[ii]) * (GG_new[ii] - GG[ii]);
      GG[ii] = in.a_mix*GG_new[ii] + (1-in.a_mix)*GG[ii];
    }
    iter_err = sqrt(iter_err);
    
    // End timing
    double toc = omp_get_wtime();
    
    // Print diagnostic
    if (verbose) {
      printf("--- iteration %d ---\n", iter_counter);
      printf("Elapsed time: %f seconds\n", toc - tic);
      printf("Residual error: %.5e\n", iter_err);
      fflush(stdout);
    }
  }
  if (verbose) printf("Done.\n");
  
  // Internal energy
  if (verbose) printf("Internal energy: %.10f\n",compute_uex(SS, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text(SS, GG, phi, SSHF, xx, in);
  write_guess(SS, GG, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
  free(bf);
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_hnc(double *GG_new, double *GG, double *SS,
		      double *bf, double *xx, input in) {

  // Static local field correction
  #pragma omp parallel for
  for (int ii=0; ii<in.nx; ii++) {

    if (xx[ii] > 0.0){

      double kmax, kmin, fu, 
	xx2, uu, uu2, ww2, ww;
   
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
	
	GG_new[ii] += (1.0 - (GG[jj] - 1.0) * (SS[jj] - 1.0) 
		       - bf[jj])*in.dx*fu/uu;
	
      }
      
      GG_new[ii] *= 3.0*in.dx/(8.0*xx[ii]);
      GG_new[ii] += bf[ii];
      
    }
    else
      GG_new[ii] = 0.0;
  }

}

void compute_bf(double *bf, double *xx, input in, bool iet){

  double scaling = 1.0;
  double ll = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double l2 = ll*ll, l3 = l2*ll, l4 = l3*ll, l5 = l4*ll, 
    l6 = l5*ll, l7 = l6*ll, l8 = l7*ll;
  double Gamma = scaling*2*l2*in.rs/in.Theta;
  double lnG = log(Gamma), lnG2 = lnG*lnG;
  double b0 = 0.258 - 0.0612*lnG + 0.0123*lnG2 - 1.0/Gamma;
  double b1 = 0.0269 + 0.0318*lnG + 0.00814*lnG2;
  double c1 = 0.498 - 0.280*lnG + 0.0294*lnG2;
  double c2 = -0.412 + 0.219*lnG - 0.0251*lnG2;
  double c3 = 0.0988 - 0.0534*lnG + 0.00682*lnG2;
  double b02 = b0*b0, b03 = b02*b0, b04 = b03*b0, b05 = b04*b0,
    b06 = b05*b0, b07 = b06*b0, b08 = b07*b0;
  double b12 = b1*b1, b13 = b12*b1, b14 = b13*b1, b15 = b14*b1,
    b16 = b15*b1, b17 = b16*b1, b18 = b17*b1;
  double b02_b12 = b02/b12, b03_b13 = b03/b13, b04_b14 = b04/b14,
    b05_b15 = b05/b15, b06_b16 = b06/b16, b07_b17 = b07/b17, 
    b08_b18 = b08/b18;
  double ff = 0.0;
  double q2,q3,q4,q5,q6,q7,q8;
  double bf1, bf2, bf3;

  if (iet){
    if (b0/b1 >= 0.0)
      ff = sqrt(M_PI)/(4.0*l2)*pow(b0/b1, 1.5);
    else{
      printf("Error: The STLS-IET scheme cannot be applied to this state point"
	     "(Gamma = %.8f) because the bridge function term diverges\n", Gamma);
      exit(EXIT_FAILURE);
    }
  }

  for (int ii=0; ii<in.nx; ii++){

    if (iet){
      q2 = xx[ii]*xx[ii];
      q3 = q2*xx[ii];
      q4 = q3*xx[ii];
      q5 = q4*xx[ii];
      q6 = q5*xx[ii];
      q7 = q6*xx[ii]; 
      q8 = q7*xx[ii];
      bf1 = -b0 + c1/16.0*(60.0*b02_b12 - 20.0*b03_b13*q2/l2 + b04_b14*q4/l4);
      bf2 = c2/64.0*(840.0*b03_b13 - 420.0*b04_b14*q2/l2 +
		     42.0*b05_b15*q4/l4 - b06_b16*q6/l6);
      bf3 = c3/256.0*(15120.0*b04_b14 - 10080.0*b05_b15*q2/l2 +
		     1512.0*b06_b16*q4/l4 - 72.0*b07_b17*q6/l6 + 
		     b08_b18*q8/l8);
      bf[ii] = scaling*ff*q2*(bf1 + bf2 + bf3)*exp(-b0*q2/(4.0*b1*l2));
    }
    else 
      bf[ii] = 0.0;
      
  }

}
