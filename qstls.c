#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "stls.h"
#include "qstls.h"


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS-HNC EQUATIONS
// -------------------------------------------------------------------

void solve_qstls(input in, bool verbose) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *psi = NULL;
  double *SS = NULL;
  double *SS_new = NULL;
  double *SSHF = NULL;

  // Solve STLS equation for initial guess
  if (verbose) printf("Solution of classical STLS for initial guess:\n");
  double a_mix_hold = in.a_mix;
  in.a_mix = 0.1;
  solve_stls(in, false, &xx, &SS, &SSHF, NULL, NULL, &phi);
  in.a_mix = a_mix_hold;
  if (verbose) printf("Done.\n");

  // Allocate arrays for qstls calculation
  alloc_qstls_arrays(in, &psi, &SS_new);  
  
  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    clock_t tic = clock();
    
    // Update auxilliary response 
    compute_psi(psi, xx, SS, in);
    
    // Update SSF
    compute_ssf(SS_new, SSHF, phi, psi, xx, in);

    // Update diagnostic
    iter_err = 0.0;
    iter_counter++;
    for (int ii=0; ii<in.nx; ii++) {
      iter_err += (SS_new[ii] - SS[ii]) * (SS_new[ii] - SS[ii]);
      SS[ii] = in.a_mix*SS_new[ii] + (1-in.a_mix)*SS[ii];
    }
    iter_err = sqrt(iter_err);
        
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
  if (verbose) printf("Internal energy: %f\n",compute_uex(SS, in));
  
  // Output to file
  /* if (verbose) printf("Writing output files...\n"); */
  /* write_text_qstls(SS, GG, xx, in); */
  /* if (verbose) printf("Done.\n"); */

  // Output to variable or free memory
  free_qstls_arrays(xx, phi, psi, SS, SS_new, SSHF);
 
 
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_qstls_arrays(input in, double **psi, double **SS_new){

  *psi = malloc( sizeof(double) * in.nx * in.nl);
  *SS_new = malloc( sizeof(double) * in.nx);

}

void free_qstls_arrays(double *xx, double *phi,
                      double *psi, double *SS,
                      double *SS_new, double *SSHF){

  free(xx);
  free(phi);
  free(psi);
  free(SSHF);
  free(SS);
  free(SS_new);

}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE AUXILLIARY RESPONSE
// -------------------------------------------------------------------

void compute_psi(double *psi, double *xx,  double *SS, input in) {

  // Temporary array to store results
  double *psil = malloc( sizeof(double) * in.nx);

  // Loop over the Matsubara frequency
  for (int ll=0; ll<in.nl; ll++){
    // Compute auxilliary response
    compute_psil(psil, xx, SS, ll, in);
    // Fill output array
    for (int ii=0; ii<in.nx; ii++){
      psi[idx2(ii,ll,in.nx)] = psil[ii];
    }
  }

  // Free memory
  free(psil);

}

void compute_psil(double *psil, double *xx, double *SS, int ll, input in) {

  double nu = (int)floor(2.0/in.dx);
  double fw, fq;

  for (int ii=0; ii<in.nx; ii++) {

    psil[ii] = 0.0;

    for (int jj=0; jj<in.nx; jj++){
      
      fw = 0.0;

      for (int kk=0; jj<in.nx; kk++){

  	fq = 0.0;

  	for (int mm=0; mm<nu; mm++){
  	  fq += psi_u(-1.0 + mm*in.dx, xx[kk], xx[jj],
  		      xx[ii], ll, in);
  	}

  	fw += fq*in.dx*psi_q(xx[kk], ll, in);
	
      }

      psil[ii] += fw*in.dx*psi_w(xx[ii], SS[ii]);

    }

    if (ll ==0) psil[ii] *= -(3.0/8.0)*in.dx;
    else psil[ii] *= -3.0/(4.0*in.Theta)*in.dx;

  }

}


/* double psi_t(double tt, double qq, double ww,  */
/* 	     double xx, int ll, input in){ */

/*   double tt2 = tt*tt, xx2 = xx*xx, qq2 = qq*qq, */
/*     txq = 2.0*xx*qq, tplT = 2*M_PI*ll*in.Theta, */
/*     tplT2 = tplT*tplT, fact = 1.0/(2*tt + ww2 - xx2); */
     
/*   if (ll == 0){ */

/*     if (tt > txq){ */
/*       return fact * ((qq2 - tt2/(4.0*xx2))* */
/* 		     log((tt + txq)/(tt - txq))  */
/* 		     + qq*tt/xx); */
/*     } */
/*     else if (tt < txq){ */
/*       return fact * ((qq2 - tt2/(4.0*xx2))* */
/* 		     log((t + txq)/(txq - tt))  */
/* 		     + qq*tt/xx); */
/*     } */
/*     else { */
/*       return fact * 2.0 * qq2; */
/*     } */

/*   } */
/*   else { */

/*     return fact * log(((txq + tt)*(txq + tt) + tplT2) */
/* 		      /((txq - tt)*(txq - tt) + tplT2)); */

/*   } */

/* } */

double psi_u(double uu, double qq, double ww,
	     double xx, int ll, input in){

  double xx2 = xx*xx, qq2 = qq*qq, ww2 = ww*ww, 
    txq = 2.0*xx*qq, tplT = 2*M_PI*ll*in.Theta, 
    tplT2 = tplT*tplT, xwu = xx*ww*uu, 
    tt = xx2 - xwu, tt2 = tt*tt, 
    fact = ww*xx/(ww2 + tt - xwu);
  
  if (ll == 0){

    if (tt > txq){
      return fact * ((qq2 - tt2/(4.0*xx2))*
		     log((tt + txq)/(tt - txq))
		     + qq*tt/xx);
    }
    else if (tt < txq){
      return fact * ((qq2 - tt2/(4.0*xx2))*
		     log((tt + txq)/(txq - tt))
		     + qq*tt/xx);
    }
    else {
      return fact * 2.0 * qq2;
    }

  }
  else {

    return fact * log(((txq + tt)*(txq + tt) + tplT2)
		      /((txq - tt)*(txq - tt) + tplT2));

  }

}

double psi_q(double qq, int ll, input in){

  double qq2 = qq*qq;

  if (ll == 0){

    return qq/(exp(qq2/in.Theta - in.mu) 
		+ exp(-qq2/in.Theta + in.mu) + 2.0);

  }
  else {

    return qq/(exp(qq2/in.Theta - in.mu) +  1.0);

  }

}

double psi_w(double ww, double SS){

  return ww * (SS - 1);
  
}

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------


// write text files with SSF and SLFC
void write_text_qstls(double *SS, double *GG, double *xx, input in){


    FILE* fid;
    
    // Output for SSF
    fid = fopen("ssf_QSTLS.dat", "w");
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
    fid = fopen("slfc_QSTLS.dat", "w");
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

