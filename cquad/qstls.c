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

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS-HNC EQUATIONS
// -------------------------------------------------------------------

void solve_qstls(input in, bool verbose) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *GG_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;
  double *psi = NULL;
  double *psi_xlw = NULL;

  // Allocate arrays
  alloc_stls_arrays(in, &xx, &phi, &GG, &GG_new, &SS, &SSHF);
  psi = malloc( sizeof(double) * in.nx * in.nl);  
  psi_xlw = malloc( sizeof(double) * in.nx * in.nl * in.nx);

  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_stls_arrays(&in, xx, phi, SSHF, verbose);
  

  /* // Initial guess for Static structure factor (SSF) and static-local field correction (SLFC) */
  /* if (strcmp(in.guess_file,"NO_FILE")==0){ */
  /*   for (int ii=0; ii < in.nx; ii++) { */
  /*     GG[ii] = 0.0; */
  /*     GG_new[ii] = 1.0; */
  /*   } */
  /*   compute_ssf(SS, SSHF, GG, phi, xx, in); */
  /* } */
  /* else { */
  /*   read_guess(SS, GG, in); */
  /* } */
   
  /* // SSF and SLFC via iterative procedure */
  /* if (verbose) printf("SSF and SLFC calculation...\n"); */
  /* double iter_err = 1.0; */
  /* int iter_counter = 0; */
  /* while (iter_counter < in.nIter && iter_err > in.err_min_iter ) { */
    
  /*   // Start timing */
  /*   double tic = omp_get_wtime(); */
    
  /*   // Update SSF */
  /*   compute_ssf(SS, SSHF, GG, phi, xx, in); */
    
  /*   // Update SLFC */
  /*   compute_slfc_hnc(GG_new, GG, SS, bf, xx, in); */
    
  /*   // Update diagnostic */
  /*   iter_err = 0.0; */
  /*   iter_counter++; */
  /*   for (int ii=0; ii<in.nx; ii++) { */
  /*     iter_err += (GG_new[ii] - GG[ii]) * (GG_new[ii] - GG[ii]); */
  /*     GG[ii] = in.a_mix*GG_new[ii] + (1-in.a_mix)*GG[ii]; */
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
  /* write_text(SS, GG, phi, SSHF, xx, in); */
  /* write_guess(SS, GG, in); */
  /* if (verbose) printf("Done.\n"); */

  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
  free(psi);
  free(psi_xlw);

}


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE ...
// -------------------------------------------------------------------

struct psi_xlw_q_params {

  double mu;
  double Theta;
  gsl_spline *t_int_sp_ptr;
  gsl_interp_accel *t_int_acc_ptr;

};

struct psi_xlw_t_params {

  double Theta;
  double ww;
  double xx;
  double ll;
  double qq;

};


void compute_psi_xlw(double *psi_xlw, double *xx, input in) {

  double err;
  size_t nevals;
  double xx2, xw, tmax, tmin;
  double *t_int  = malloc( sizeof(double) * in.nx);

  // Declare accelerator and spline objects
  gsl_spline *t_int_sp_ptr;
  gsl_interp_accel *t_int_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  t_int_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  t_int_acc_ptr = gsl_interp_accel_alloc();

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
  = gsl_integration_cquad_workspace_alloc(100);

  // Loop over xx (wave-vector)
  for (int ii=0; ii<in.nx; ii++){

    // Loop over ll (Matsubara frequencies)
    for (int ll=0; ll<in.nl; ll++){
      
      // Integration function
      gsl_function ft_int, fq_int;
      if (ll == 0){
	ft_int.function = &psi_x0w_t;
	fq_int.function = &psi_x0w_q;
      }
      else {
	ft_int.function = &psi_xlw_t;
	fq_int.function = &psi_xlw_q;
      }

      // Loop over w (wave-vector)
      for (int jj=0; jj<in.nx; jj++){

	// Integration limits for the integration over t 
	xx2 = xx[ii]*xx[ii];
	xw = xx[ii]*xx[jj];
	tmin = xx2 - xw + in.dx; // +in.dx was added to avoid overflows
	tmax = xx2 + xw;

	// Construct integrand for the integral over q
	for (int kk=0; kk<in.nx; ii++) {

	  if (xx[kk] > 0.0){
	    
	    // Integration over t
	    struct psi_xlw_t_params ppt = {in.Theta,xx[jj],xx[ii],ll,xx[kk]};
	    ft_int.params = &ppt;
	    gsl_integration_cquad(&ft_int,
				  tmin, tmax,
				  0.0, 1e-5,
				  wsp,
				  &t_int[kk], &err, &nevals);
	  }

	  else t_int[kk] = 0.0;

	}
	gsl_spline_init(t_int_sp_ptr, xx, t_int, in.nx);
      
	// Integral over q
	struct psi_xlw_q_params ppq = {in.mu,in.Theta,t_int_sp_ptr,t_int_acc_ptr};
	fq_int.params = &ppq;
	gsl_integration_cquad(&fq_int,
			      xx[0], xx[in.nx-1],
			      0.0, 1e-5,
			      wsp,
			      &psi_xlw[idx3(ii, ll, jj, in.nx, in.nl)],
			      &err, &nevals);
      }
    }
  }


  // Free memory
  free(t_int);
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(t_int_sp_ptr);
  gsl_interp_accel_free(t_int_acc_ptr);

}

double psi_x0w_t(double tt, void* pp) {

  struct psi_xlw_t_params* params = (struct psi_xlw_t_params*)pp;
  double xx = (params->xx);
  double qq = (params->qq);
  double ww = (params->ww);
  double xx2 = xx*xx, ww2 = ww*ww, qq2 = qq*qq, 
    tt2 = tt*tt, txq = 2.0*xx*qq;
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

double psi_x0w_q(double qq, void* pp) {

  struct psi_xlw_q_params* params = (struct psi_xlw_q_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  gsl_spline* t_int_sp_ptr = (params->t_int_sp_ptr);
  gsl_interp_accel* t_int_acc_ptr = (params->t_int_acc_ptr);
  double qq2 = qq*qq;
  double fft = gsl_spline_eval(t_int_sp_ptr, qq, t_int_acc_ptr);

  return qq/(exp(qq2/Theta - mu) + exp(-qq2/Theta + mu) + 2.0)*fft;

}

double psi_xlw_t(double tt, void* pp) {
  
  struct psi_xlw_t_params* params = (struct psi_xlw_t_params*)pp;
  double xx = (params->xx);
  double qq = (params->qq);
  double ww = (params->ww);
  double ll = (params->ll);
  double Theta = (params->Theta);
  double xx2 = xx*xx, ww2 = ww*ww, txq = 2.0*xx*qq, 
    tplT = 2.0*M_PI*ll*Theta, tplT2 = tplT*tplT, 
    txqpt = txq + tt, txqmt = txq - tt,
    txqpt2 = txqpt*txqpt, txqmt2 = txqmt*txqmt,
    logarg = (txqpt2 + tplT2)/(txqmt2 + tplT2);

  if (xx == 0 || qq == 0)
    return 0;
  else
    return 1.0/(2.0*tt + ww2 - xx2)*log(logarg);

}

double psi_xlw_q(double qq, void* pp) {
  
  struct psi_xlw_q_params* params = (struct psi_xlw_q_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  gsl_spline* t_int_sp_ptr = (params->t_int_sp_ptr);
  gsl_interp_accel* t_int_acc_ptr = (params->t_int_acc_ptr);
  double qq2 = qq*qq;
  double fft = gsl_spline_eval(t_int_sp_ptr, qq, t_int_acc_ptr);

  return qq/(exp(qq2/Theta - mu) + 1.0)*fft;

}


// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A THREE-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx3(int xx, int yy, int zz,
         int x_size, int y_size) {
  return (zz * x_size * y_size) + (yy * x_size) + xx;
}

