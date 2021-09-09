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
#include "stls_iet.h"

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS-IET EQUATIONS
// -------------------------------------------------------------------

void solve_stls_iet(input in, bool verbose) {

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
  compute_bf(bf, xx, in);

  // Initial guess for Static structure factor (SSF) and static-local field correction (SLFC)
  if (strcmp(in.guess_file,"NO_FILE")==0){
    for (int ii=0; ii < in.nx; ii++) {
      GG[ii] = 0.0;
      GG_new[ii] = 1.0;
    }
    compute_ssf_static(SS, SSHF, GG, phi, xx, in);
  }
  else {
    read_guess_static(SS, GG, in);
  }
   
  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    double tic = omp_get_wtime();
    
    // Update SSF
    compute_ssf_static(SS, SSHF, GG, phi, xx, in);
    
    // Update SLFC
    compute_slfc_iet(GG_new, GG, SS, bf, xx, in);
    
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
  if (verbose) printf("Internal energy: %.10f\n",compute_uex(SS, xx, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text_static(SS, GG, phi, SSHF, xx, in);
  write_bf_static(bf, xx, in);
  write_guess_static(SS, GG, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
  free(bf);
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

struct slfcu_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *slfc_sp_ptr;
  gsl_interp_accel *slfc_acc_ptr;
  gsl_spline *GGu_sp_ptr;
  gsl_interp_accel *GGu_acc_ptr;
  gsl_spline *bf_sp_ptr;
  gsl_interp_accel *bf_acc_ptr;

};

struct slfcw_params {

  double xx;
  double uu;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};


void compute_slfc_iet(double *GG_new, double *GG, double *SS,
                      double *bf, double *xx, input in) {

  double err;
  size_t nevals;
  double wmax, wmin, GG_tmp;
  double *GGu  = malloc( sizeof(double) * in.nx);

  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  gsl_spline *slfc_sp_ptr;
  gsl_interp_accel *slfc_acc_ptr;
  gsl_spline *GGu_sp_ptr;
  gsl_interp_accel *GGu_acc_ptr;
  gsl_spline *bf_sp_ptr;
  gsl_interp_accel *bf_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  slfc_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  slfc_acc_ptr = gsl_interp_accel_alloc();
  GGu_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  GGu_acc_ptr = gsl_interp_accel_alloc();
  bf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  bf_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);    
  gsl_spline_init(slfc_sp_ptr, xx, GG, in.nx);    
  gsl_spline_init(bf_sp_ptr, xx, bf, in.nx);    

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
  = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function fu_int, fw_int;
  fu_int.function = &slfc_u;
  fw_int.function = &slfc_w;

  // STLS component of the static local field correction
  compute_slfc(GG_new, SS, xx, in);

  // Non-STLS component of the static local field correction
  // Integration over u
  for (int ii=0; ii<in.nx; ii++) {

    if (xx[ii] > 0.0){

      // Integration over w
      for (int jj=0; jj<in.nx; jj++){
	
	if (xx[jj] >  0.0) {

	  struct slfcw_params slfcwp = {xx[ii], xx[jj], 
					ssf_sp_ptr, ssf_acc_ptr};
	  wmin = xx[jj] - xx[ii];
	  if (wmin < 0.0) wmin = -wmin;
	  // NOTE: The upper cutoff is at qm - dq for numerical reasons;
	  // The quadrature formula attemps a tiny extrapolation which causes 
	  // the interpolation routine to crash. 
	  wmax = GSL_MIN(xx[in.nx-2], xx[ii]+xx[jj]);
	  fw_int.params = &slfcwp;
	  gsl_integration_cquad(&fw_int,
				wmin, wmax,
				0.0, 1e-5,
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
				    bf_sp_ptr, bf_acc_ptr};
      fu_int.params = &slfcup;
      gsl_integration_cquad(&fu_int,
			    xx[0], xx[in.nx-1],
			    0.0, 1e-5,
			    wsp,
			    &GG_tmp, &err, &nevals);
      
      GG_new[ii] += 3.0/(8.0*xx[ii])*GG_tmp + bf[ii];

    }
    else 
      GG_new[ii] = 0.0;

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
  gsl_spline_free(bf_sp_ptr);
  gsl_interp_accel_free(bf_acc_ptr);

}

double slfc_u(double uu, void* pp) {

  struct slfcu_params* params = (struct slfcu_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  gsl_spline* slfc_sp_ptr = (params->slfc_sp_ptr);
  gsl_interp_accel* slfc_acc_ptr = (params->slfc_acc_ptr);
  gsl_spline* GGu_sp_ptr = (params->GGu_sp_ptr);
  gsl_interp_accel* GGu_acc_ptr = (params->GGu_acc_ptr);
  gsl_spline* bf_sp_ptr = (params->bf_sp_ptr);
  gsl_interp_accel* bf_acc_ptr = (params->bf_acc_ptr);

  if (uu > 0.0)
    return (1.0/uu) * gsl_spline_eval(GGu_sp_ptr, uu, GGu_acc_ptr)
      *(-gsl_spline_eval(bf_sp_ptr, uu, bf_acc_ptr) 
	- (gsl_spline_eval(ssf_sp_ptr, uu, ssf_acc_ptr) - 1.0)
	*(gsl_spline_eval(slfc_sp_ptr, uu, slfc_acc_ptr) - 1.0));
  else 
    return 0;
}

double slfc_w(double ww, void* pp) {

  struct slfcw_params* params = (struct slfcw_params*)pp;
  double xx = (params->xx);
  double uu = (params->uu);
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  double ww2 = ww*ww, xx2 = xx*xx, uu2 = uu*uu;
    
  return (ww2 - uu2 - xx2)*ww
    *(gsl_spline_eval(ssf_sp_ptr, ww, ssf_acc_ptr) - 1.0);

}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE BRIDGE FUNCTION TERM
// -------------------------------------------------------------------

void compute_bf(double *bf, double *xx, input in){

  if (in.theory_id == 2 || in.theory_id == 7)
    bf_hnc(bf, xx, in);
  else if (in.theory_id == 3 || in.theory_id == 8)
    bf_ocp_ioi(bf, xx, in);
  else if (in.theory_id == 4 || in.theory_id == 9)
    bf_ocp_lct(bf, xx, in);
  else if (in.theory_id == 5)
    bf_rescaled_ocp_lct(bf, xx, in);
  else{
    printf("Error: unknown theory to be compute the bridge function."
           "Choose between: STLS-IET-HNC, STLS-IET-IOI, STLS-IET-LCT, STLS-RIET-LCT,"
	   "QSTLS-IET-HNC, QSTLS-IET-IOI, QSTLS-IET-LCT\n");
    exit(EXIT_FAILURE);
  }

}


// HNC bridge function  
void bf_hnc(double *bf, double *xx, input in){

  for (int ii=0; ii<in.nx; ii++){
      bf[ii] = 0.0;      
  }

}
 
// Bridge function from the parameterization of Ichimaru and collaborators
void bf_ocp_ioi(double *bf, double *xx, input in){

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

  if (b0/b1 <= 0.0 || Gamma < 5.25 || Gamma > 171.8){

    printf("Error: The STLS-IET scheme cannot be applied to this state point"
	   " because Gamma = %.8f falls outside the range of validty of the"
	   " bridge function parameterization\n", Gamma);
    exit(EXIT_FAILURE);

  }

  ff = sqrt(M_PI)/(4.0*l2)*pow(b0/b1, 1.5);
  for (int ii=0; ii<in.nx; ii++){

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
 
}

// Bridge function from the parameterization of Lucco Castello and Tolias
struct bfr_params {

  double Gamma;

};


void bf_ocp_lct(double *bf, double *xx, input in){

  double ll = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double l2 = ll*ll;
  double Gamma = 2*l2*in.rs/in.Theta;

  double err;
  
  // Integration workspace
  gsl_integration_workspace *wsp 
    = gsl_integration_workspace_alloc(1000);
  gsl_integration_workspace *wspc 
    = gsl_integration_workspace_alloc(1000);
  gsl_integration_qawo_table *qtab
    = gsl_integration_qawo_table_alloc(0.0, 1.0, GSL_INTEG_SINE, 1000);
  
  // Integration function
  gsl_function ff_int;
  struct bfr_params rbfrp = {Gamma};
  ff_int.function = &rbfr;
  ff_int.params = &rbfrp;
  
  // Bridge function term (B(q)/U(q))
  for (int ii = 0; ii < in.nx; ii++) {

    // Set wave-vector (divide xx[ii] by ll to convert to Wigner-Seitz units)
    gsl_integration_qawo_table_set(qtab, xx[ii]/ll, 1.0, GSL_INTEG_SINE);

    // Fourier transform
    gsl_integration_qawf(&ff_int,
    			 0.0,
    			 1e-10, 1000,
    			 wsp, wspc,
    			 qtab,
    			 &bf[ii], &err);
    bf[ii] *= xx[ii]/Gamma/ll;

  }

  // Free memory
  gsl_integration_workspace_free(wsp);
  gsl_integration_workspace_free(wspc); 
  gsl_integration_qawo_table_free(qtab);

 

}


void bf_rescaled_ocp_lct(double *bf, double *xx, input in){

  double ll = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double l2 = ll*ll;
  double Gamma = 2*l2*in.rs/sqrt(1 + in.Theta*in.Theta);

  double err;
  
  // Integration workspace
  gsl_integration_workspace *wsp 
    = gsl_integration_workspace_alloc(1000);
  gsl_integration_workspace *wspc 
    = gsl_integration_workspace_alloc(1000);
  gsl_integration_qawo_table *qtab
    = gsl_integration_qawo_table_alloc(0.0, 1.0, GSL_INTEG_SINE, 1000);
  
  // Integration function
  gsl_function ff_int;
  struct bfr_params rbfrp = {Gamma};
  ff_int.function = &rbfr;
  ff_int.params = &rbfrp;
  
  // Bridge function term (B(q)/U(q))
  for (int ii = 0; ii < in.nx; ii++) {

    // Set wave-vector (divide xx[ii] by ll to convert to Wigner-Seitz units)
    gsl_integration_qawo_table_set(qtab, xx[ii]/ll, 1.0, GSL_INTEG_SINE);

    // Fourier transform
    gsl_integration_qawf(&ff_int,
    			 0.0,
    			 1e-10, 1000,
    			 wsp, wspc,
    			 qtab,
    			 &bf[ii], &err);
    bf[ii] *= xx[ii]/Gamma/ll;

  }

  // Free memory
  gsl_integration_workspace_free(wsp);
  gsl_integration_workspace_free(wspc); 
  gsl_integration_qawo_table_free(qtab);

 

}

double  rbfr(double rr, void *pp){

  struct bfr_params *params = (struct bfr_params*)pp;
  double Gamma = (params->Gamma);
  double Gamma1_6 = pow(Gamma, 1./6.);
  double lnG = log(Gamma), lnG2 = lnG*lnG; 
  double  lnG3 = lnG2*lnG, lnG4 = lnG3*lnG;
  double a0 = Gamma * (0.076912 - 0.10465*lnG + 0.0056629*lnG2
  		       + 0.00025656*lnG3);

  double a2 = Gamma * (0.068045 - 0.036952*lnG + 0.048818*lnG2
  		       - 0.0048985*lnG3);

  double a3 = Gamma * (-0.30231 + 0.30457*lnG - 0.11424*lnG2
  		       + 0.0095993*lnG3);

  double a4 = Gamma * (0.25111 - 0.26800*lnG + 0.082268*lnG2
  		       - 0.0064960*lnG3);

  double a5 = Gamma * (-0.061894 + 0.066811*lnG - 0.019140*lnG2
  		       + 0.0014743*lnG3);

  double c0 = Gamma * (0.25264 - 0.31615*lnG + 0.13135*lnG2 
		       - 0.023044*lnG3 + 0.0014666*lnG4);
 
  double c1 = Gamma1_6 * (-12.665 + 20.802*lnG - 9.6296*lnG2 
			  + 1.7889*lnG3 - 0.11810*lnG4); 

  double c2 = Gamma1_6 * (15.285 - 14.076*lnG + 5.7558*lnG2 
			  - 1.0188*lnG3 + 0.06551*lnG4); 

  double c3 = Gamma1_6 * (35.330 - 40.727*lnG + 16.690*lnG2 
			  - 2.8905*lnG3 + 0.18243*lnG4);

  double r2, r3, r4, r5, rshift;
  double bsr, blr, ff;

  if (Gamma < 5.0){

    printf("Error: The STLS-IET scheme cannot be applied to this state point"
	   " because for Gamma = %.8f the bridge function parameterization"
	   " is not applicable (Gamma must be larger or equal than 5.0)\n", Gamma);
    exit(EXIT_FAILURE);

  }
     
  r2 = rr*rr;
  r3 = r2*rr;
  r4 = r3*rr;
  r5 = r4*rr;
  rshift = rr - 1.44;

  // Short range fit (cut to avoid numerical problems)
  bsr = a0 + a2*r2 + a3*r3 + a4*r4 + a5*r5;

  // Long range fit
  blr = c0 * exp(-c1*rshift) * exp(-0.3*r2) * ( cos(c2*rshift) + c3*exp(-3.5*rshift) );

  // Full range fit
  ff = 0.5 * ( 1.0 + erf(5.0*(rr - 1.50)) );

  return rr*((1 - ff)*bsr + ff*blr);

}

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------


// write text files for output
void write_bf_static(double *bf, double *xx, input in){


  FILE* fid;

  // Output for the bridge function
  char out_name[100];
  sprintf(out_name, "bf_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    perror("Error while creating the output file for the bridge function");
    exit(EXIT_FAILURE);
  }
  for (int ii = 0; ii < in.nx; ii++)
    fprintf(fid, "%.8e %.8e\n", xx[ii], bf[ii]);

  fclose(fid);

}
