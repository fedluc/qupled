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
  double *SS_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;
  double *psi = NULL;
  double *psi_xlw = NULL;

  // Allocate arrays
  alloc_stls_arrays(in, &xx, &phi, &GG, &SS_new, &SS, &SSHF);
  psi = malloc( sizeof(double) * in.nx * in.nl);  
  psi_xlw = malloc( sizeof(double) * in.nx * in.nl * in.nx);

  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_stls_arrays(&in, xx, phi, SSHF, verbose);
  if (verbose) printf("Fixed component of the auxiliary response function: ");  
  compute_psi_xlw(psi_xlw, xx, in);
  if (verbose) printf("Done.\n");

  // Initial guess for Static structure factor (SSF) and for auxilliary response (DLFC)
  if (strcmp(in.guess_file,"NO_FILE")==0){
    for (int ii=0; ii<in.nx; ii++){
      for (int ll=0; ll<in.nl; ll++){
	psi[idx2(ii,ll,in.nx)] = 0.0;
      }
    }
    compute_ssf_qstls(SS, SSHF, psi, phi, xx, in);
  }
  else {
    read_guess(SS, GG, in);
  }
   
  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    double tic = omp_get_wtime();
    
    // Update auxiliary function
    compute_psi(psi, psi_xlw, SS, xx, in);
    
    // Update SSF
    compute_ssf_qstls(SS_new, SSHF, psi, phi, xx, in);
    
    // Update diagnostic
    iter_err = 0.0;
    iter_counter++;
    for (int ii=0; ii<in.nx; ii++) {
      iter_err += (SS_new[ii] - SS[ii]) * (SS_new[ii] - SS[ii]);
      GG[ii] = in.a_mix*SS_new[ii] + (1-in.a_mix)*SS[ii];
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
  if (verbose) printf("Internal energy: %.10f\n",compute_uex(SS, xx, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text(SS, GG, phi, SSHF, xx, in);
  write_guess(SS, GG, in);
  // THE OUTPUT ROUTINES MUST BE MODIFIED
  if (verbose) printf("Done.\n");

  // Free memory
  free_stls_arrays(xx, phi, GG, SS_new, SS, SSHF);
  free(psi);
  free(psi_xlw);

}


// ------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FIXED COMPONENT OF THE AUXILIARY RESPONSE
// ------------------------------------------------------------------------

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

  // Loop over xx (wave-vector)
  #pragma omp parallel for
  for (int ii=0; ii<in.nx; ii++){    

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
	tmax = xx2 + xw + in.dx;

	// Construct integrand for the integral over q
	for (int kk=0; kk<in.nx; kk++) {

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

    // Free memory
    free(t_int);
    gsl_integration_cquad_workspace_free(wsp);
    gsl_spline_free(t_int_sp_ptr);
    gsl_interp_accel_free(t_int_acc_ptr);

  }


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

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE CHANGING COMPONENT OF THE AUXILIARY RESPONSE
// ---------------------------------------------------------------------------

struct psiw_params {

  gsl_spline *qt_int_sp_ptr;
  gsl_interp_accel *qt_int_acc_ptr;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};


void compute_psi(double *psi, double *psi_xlw, double *SS, 
		 double *xx, input in){

  double err;
  size_t nevals;
  double *qt_int  = malloc( sizeof(double) * in.nx);
  
  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *qt_int_sp_ptr;
  gsl_interp_accel *qt_int_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  qt_int_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  qt_int_acc_ptr = gsl_interp_accel_alloc();
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  
  // Interpolate SSF
  gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);
	
  // Integration function
  gsl_function psiw_int;
  psiw_int.function = &psiw;
  
  for (int ii=0; ii<in.nx; ii++){
    for (int ll=0; ll<in.nl; ll++){

      // Interpolate solution of q-t integration
      for (int jj=0; jj<in.nl; jj++){
	qt_int[jj] = psi_xlw[idx3(ii,ll,jj,in.nx,in.nl)];
      }
      gsl_spline_init(qt_int_sp_ptr, xx, qt_int, in.nx);

      // Integral over w 
      struct psiw_params ppw = {qt_int_sp_ptr,qt_int_acc_ptr,
				ssf_sp_ptr, ssf_acc_ptr};
      psiw_int.params = &ppw;
      gsl_integration_cquad(&psiw_int,
			    xx[0], xx[in.nx-1],
			    0.0, 1e-5,
			    wsp,
			    &psi[idx2(ii,ll,in.nx)], 
			    &err, &nevals);

    }
  }  

}


double psiw(double ww, void* pp) {
  
  struct psiw_params* params = (struct psiw_params*)pp;
  gsl_spline* qt_int_sp_ptr = (params->qt_int_sp_ptr);
  gsl_interp_accel* qt_int_acc_ptr = (params->qt_int_acc_ptr);
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);

  return ww*gsl_spline_eval(qt_int_sp_ptr, ww, qt_int_acc_ptr)
    *(gsl_spline_eval(ssf_sp_ptr, ww, ssf_acc_ptr) - 1.0);

}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_qstls(double *SS, double *SSHF, double *psi,
		       double *phi, double *xx, input in){

  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double ff = 4*lambda*in.rs/M_PI;
  double xx2, BB, BB_tmp, BB_den, psixl, phixl;
  
  for (int ii=0; ii<in.nx; ii++){
    
    if (xx[ii] > 0.0){
      
      xx2 = xx[ii]*xx[ii];
      BB = 0.0;
      
      for (int ll=0; ll<in.nl; ll++){
	
	psixl = psi[idx2(ii,ll,in.nx)];
	phixl = phi[idx2(ii,ll,in.nx)];
	BB_den = 1.0 + ff/xx2*(1.0 - psixl/phixl)*phixl;
	BB_tmp = phixl*phixl*(1.0 - psixl/phixl)/BB_den;
	if (ll>0) BB_tmp *= 2.0;
	BB += BB_tmp;
	
      }
      
      SS[ii] = SSHF[ii] - 3.0/2.0*ff/xx2*in.Theta*BB;
      
    }
    else
      SS[ii] = 0.0;
    
  }
  
}



// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A THREE-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx3(int xx, int yy, int zz,
         int x_size, int y_size) {
  return (zz * x_size * y_size) + (yy * x_size) + xx;
}

