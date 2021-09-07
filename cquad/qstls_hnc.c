#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <string.h>
#include "solvers.h"
#include "stls.h"
#include "qstls.h"
#include "qstls_hnc.h"

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE QSTLS EQUATIONS
// -------------------------------------------------------------------

void solve_qstls_hnc(input in, bool verbose) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *SS_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;

  // Arrays for QSTLS solution
  double *psi = NULL;
  bool psi_xluw_init = true;

  // Allocate arrays
  // Note: GG is not needed for QSTLS, but we keep it here so that
  // we can reuse some stls routines 
  alloc_stls_arrays(in, &xx, &phi, &GG, &SS_new, &SS, &SSHF);
  psi = malloc( sizeof(double) * in.nx * in.nl);  

  // Initialize STLS arrays that are not modified by the iterative procedure
  init_fixed_stls_arrays(&in, xx, phi, SSHF, verbose);

  // Initial guess
  /* if (strcmp(in.guess_file,"NO_FILE")==0){ */
  /*   for (int ii=0; ii<in.nx; ii++){ */
  /*     for (int ll=0; ll<in.nl; ll++){ */
  /* 	psi[idx2(ii,ll,in.nx)] = 0.0; */
  /*     } */
  /*   } */
  /*   compute_ssf_dynamic(SS, SSHF, psi, phi, xx, in); */
  /* } */
  /* else { */
  /*   read_guess_dynamic(SS, psi_xlw, in, &psi_xlw_init); */
  /* } */

  // Initialize QSTLS arrays that are not modified by the iterative procedure
  if (psi_xluw_init){
    if (verbose) printf("Fixed component of the auxiliary response function: ");
    compute_psi_xluw(xx, in);
    if (verbose) printf("Done.\n");
  }
 
  /* // SSF and SLFC via iterative procedure */
  /* if (verbose) printf("SSF calculation...\n"); */
  /* double iter_err = 1.0; */
  /* int iter_counter = 0; */
  /* while (iter_counter < in.nIter && iter_err > in.err_min_iter ) { */
    
  /*   // Start timing */
  /*   double tic = omp_get_wtime(); */
    
  /*   // Update auxiliary function */
  /*   compute_psi_hnc(psi, psi_xlw, SS, xx, in); */
    
  /*   // Update SSF */
  /*   compute_ssf_dynamic(SS_new, SSHF, psi, phi, xx, in); */
    
  /*   // Update diagnostic */
  /*   iter_err = 0.0; */
  /*   iter_counter++; */
  /*   for (int ii=0; ii<in.nx; ii++) { */
  /*     iter_err += (SS_new[ii] - SS[ii]) * (SS_new[ii] - SS[ii]); */
  /*     SS[ii] = in.a_mix*SS_new[ii] + (1-in.a_mix)*SS[ii]; */
  /*   } */
  /*   iter_err = sqrt(iter_err); */
    
  /*   // End timing */
  /*   double toc = omp_get_wtime(); */
    
  /*   // Print diagnostic */
  /*   if (verbose) { */
  /*     printf("--- iteration %d ---\n", iter_counter); */
  /*     printf("Elapsed time: %f seconds\n", toc - tic); */
  /*     printf("Residual error: %.5e\n", iter_err); */
  /*     fflush(stdout); */
  /*   } */
  /* } */
  /* if (verbose) printf("Done.\n"); */
  
  /* // Internal energy */
  /* if (verbose) printf("Internal energy: %.10f\n",compute_uex(SS, xx, in)); */
  
  /* // Output to file */
  /* if (verbose) printf("Writing output files...\n"); */
  /* write_text_dynamic(SS, psi, phi, SSHF, xx, in); */
  /* write_guess_dynamic(SS, psi_xlw, in); */
  /* if (verbose) printf("Done.\n"); */

  // Free memory
  free_stls_arrays(xx, phi, GG, SS_new, SS, SSHF);
  free(psi);

}


// ------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FIXED COMPONENT OF THE AUXILIARY RESPONSE
// ------------------------------------------------------------------------

struct psi_xluw_params {

  double mu;
  double Theta;
  double ll;
  double xx;
  double uu;
  double ww;

};


void compute_psi_xluw(double *xx, input in) {

  // Parallel calculations
  #pragma omp parallel
  {  

    double err;
    size_t nevals;
    double psi_xluw;
     
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);

    // Loop over xx (wave-vector)
    #pragma omp for // Distribute loops over the threads
    for (int ii=0; ii<in.nx; ii++){

      // Open binary file for output
      char out_name[100];
      sprintf(out_name, "psi_fixed_theta%.3f_xx%.5f_%s.bin", in.Theta, xx[ii], in.theory);
      FILE *fid = NULL;
      fid = fopen(out_name, "wb");
      if (fid == NULL) {
	fprintf(stderr,"Error while creating file for fixed component of the auxilliary response function");
	exit(EXIT_FAILURE);
      }

      
      // Loop over ll (Matsubara frequencies)
      for (int ll=0; ll<in.nl; ll++){
	
    	// Integration function
    	gsl_function ff_int;
    	if (ll == 0){
    	  ff_int.function = &psi_x0uw_y;
    	}
    	else {
    	  ff_int.function = &psi_xluw_y;
    	}
	
    	// Loop over u (wave-vector)
    	for (int jj=0; jj<in.nx; jj++){

  	  // Loop over w
  	  for (int kk=0; kk<in.nx; kk++) {
	    
  	    // Integral over y
  	    if (xx[ii] == 0.0 || xx[jj] == 0.0 || xx[kk] == 0.0){
  	      // No need to compute the integral over y in this cases
  	      psi_xluw = 0.0;
  	    }
  	    else{
  	      // Compute the integral and store it
  	      struct psi_xluw_params pp = {in.mu,in.Theta,ll,xx[ii],xx[jj],xx[kk]};
  	      ff_int.params = &pp;
  	      gsl_integration_cquad(&ff_int,
  				    xx[0], xx[in.nx-1],
  				    0.0, 1e-5,
  				    wsp,
  				    &psi_xluw,
  				    &err, &nevals);
           }
	    
	    // Write result to output file
	    fwrite(&psi_xluw, sizeof(double), 1, fid);
	    
	  }
    	}
      }

      // Close binary file for output
      fclose(fid);

    }

    // Free memory
    gsl_integration_cquad_workspace_free(wsp);
 
  }
  
}

double psi_x0uw_y(double yy, void* pp) {

  struct psi_xluw_params* params = (struct psi_xluw_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  double xx = (params->xx);
  double uu = (params->uu);
  double ww = (params->ww);
  double xx2 = xx*xx, ww2 =  ww*ww, uu2 = uu*uu, fxy = 4.0*xx*yy, 
    yy2 = yy*yy, umwpx = uu2 - ww2 + xx2,
    logarg = (umwpx + fxy)/(umwpx - fxy);

  if (logarg < 0.0) logarg = -logarg;

  if (xx == 0 || uu == 0 || ww == 0)
    return 0;
  else
    return yy/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0)*
      ((yy2 - umwpx*umwpx/(16.0*xx2))*log(logarg) + (yy/xx)*umwpx/2.0);

}


double psi_xluw_y(double yy, void* pp) {
  
  struct psi_xluw_params* params = (struct psi_xluw_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  double ll = (params->ll);
  double xx = (params->xx);
  double uu = (params->uu);
  double ww = (params->ww);
  double xx2 = xx*xx, ww2 =  ww*ww, uu2 = uu*uu, fxy = 4.0*xx*yy, 
    yy2 = yy*yy, fplT = 4.0*M_PI*ll*Theta, fplT2 = fplT*fplT,
    umwpx = uu2 - ww2 + xx2,
    logarg = ((umwpx + fxy)*(umwpx + fxy) + fplT2)/
             ((umwpx - fxy)*(umwpx - fxy) + fplT2);

  if (xx == 0 || uu == 0 || ww == 0)
    return 0;
  else
    return yy/(exp(yy2/Theta - mu) + 1.0)*log(logarg);

}

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE CHANGING COMPONENT OF THE AUXILIARY RESPONSE
// ---------------------------------------------------------------------------

struct psiu_params {

  gsl_spline *uint_sp_ptr;
  gsl_interp_accel *uint_acc_ptr;

};


struct psiw_params {

  gsl_spline *wint_sp_ptr;
  gsl_interp_accel *wint_acc_ptr;

};

void compute_psi_hnc(double *psi_new, double *psi, double *phi, 
		     double *SS, double *xx, input in){

  double err;
  size_t nevals;
  double wmax, wmin;
  double *wint  = malloc( sizeof(double) * in.nx);
  double *uint  = malloc( sizeof(double) * in.nx);

  // Declare accelerator and spline objects
  gsl_spline *wint_sp_ptr;
  gsl_interp_accel *wint_acc_ptr;
  gsl_spline *uint_sp_ptr;
  gsl_interp_accel *uint_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  wint_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  wint_acc_ptr = gsl_interp_accel_alloc();
  uint_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  uint_acc_ptr = gsl_interp_accel_alloc();
   
  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
	
  // Integration function
  gsl_function fu_int, fwint;
  fu_int.function = &psi_u_hnc;
  fwint.function = &psi_w_hnc;

  // Loop over xx (wave-vector)
  for (int ii=0; ii<in.nx; ii++){

    // Open binary file with the fixed component of the auxilliary response function
    char out_name[100];
    sprintf(out_name, "psi_fixed_theta%.3f_xx%.5f_%s.bin", in.Theta, xx[ii], in.theory);
    FILE *fid = NULL;
    fid = fopen(out_name, "rb");
    if (fid == NULL) {
      fprintf(stderr,"Error while reading file for fixed component of the auxilliary response function");
      exit(EXIT_FAILURE);
    }

    // Loop over the Matsubara frequencies
    for (int ll=0; ll<in.nl; ll++){

      // Loop over u (wave-vector)
      for (int jj=0; jj<in.nx; jj++){

  	// Construct integrand over w
  	for (int kk=0; kk<in.nx; kk++) {
  	  fread(&wint[kk], sizeof(double), 1, fid);
	  wint[kk] *= xx[kk]*SS[kk];
  	}
	gsl_spline_init(wint_sp_ptr, xx, wint, in.nx);

	// Compute integral over w
	struct psiw_params psiwp =  {wint_sp_ptr, wint_acc_ptr};
  	wmin = xx[jj] - xx[ii];
  	if (wmin < 0.0) wmin = -wmin;
  	// NOTE: The upper cutoff is at qm - dq for numerical reasons;
  	// The quadrature formula attemps a tiny extrapolation which causes
  	// the interpolation routine to crash.
  	wmax = GSL_MIN(xx[in.nx-2], xx[ii]+xx[jj]);
	fwint.params = &psiwp;
	gsl_integration_cquad(&fwint,
			      wmin, wmax,
			      0.0, 1e-5,
			      wsp,
			      &uint[jj], &err, &nevals);

	// Construct integrand over u
	if (xx[jj] > 0.0){
	  uint[jj] *= (1.0/xx[jj])
	    *(1 - (psi[idx2(jj,ll,in.nx)]/phi[idx2(jj,ll,in.nx)])
	     *(SS[jj]-1));
	} 
	else{
	  uint[jj] = 0.0;
	}

      }

      // Interpolate integrand over u
      gsl_spline_init(uint_sp_ptr, xx, uint, in.nx);    

      // Compute integral over u
      struct psiu_params psiup = {uint_sp_ptr, uint_acc_ptr};
      fu_int.params = &psiup;
      gsl_integration_cquad(&fu_int,
			    xx[0], xx[in.nx-1],
			    0.0, 1e-5,
			    wsp,
			    &psi_new[idx2(ii,ll,in.nx)], &err, &nevals);

    }
  }

  // Free memory
  free(wint);
  free(uint);
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(uint_sp_ptr);
  gsl_interp_accel_free(uint_acc_ptr);
  gsl_spline_free(wint_sp_ptr);
  gsl_interp_accel_free(wint_acc_ptr);

}


double psi_u_hnc(double uu, void* pp) {

  struct psiu_params* params = (struct psiu_params*)pp;
  gsl_spline* uint_sp_ptr = (params->uint_sp_ptr);
  gsl_interp_accel* uint_acc_ptr = (params->uint_acc_ptr);
 
  return gsl_spline_eval(uint_sp_ptr, uu, uint_acc_ptr);

}

double psi_w_hnc(double ww, void* pp) {

  struct psiw_params* params = (struct psiw_params*)pp;
  gsl_spline* wint_sp_ptr = (params->wint_sp_ptr);
  gsl_interp_accel* wint_acc_ptr = (params->wint_acc_ptr);
 
  return gsl_spline_eval(wint_sp_ptr, ww, wint_acc_ptr);
}


