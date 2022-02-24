#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include "solvers.h"
#include "utils.h"
#include "chemical_potential.h"
#include "stls.h"
#include "vs_stls.h"

// -------------------------------------------------------------------
// FUNCTION USED TO SOLVE THE VS-STLS SCHEME (AT ZERO TEMPERATURE)
// -------------------------------------------------------------------

void solve_vs_stls(input in, bool verbose) {

  // VS parameters
  double alpha = 0.5;
  double nrs = 10;
  double rs = in.rs;
  double drs;
  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *GG_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;

  // Arrays for VS-STLS solution
  double *uint = NULL;
  double *rsArray = NULL;
  
  // Allocate arrays
  alloc_stls_arrays(in, &xx, &phi, &GG, &GG_new, &SS, &SSHF);
  alloc_vs_stls_arrays(in, nrs, &uint, &rsArray);
  
  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_vs_stls_arrays(&in, nrs, xx, rsArray, verbose);
  
  // Initial guess (not sure how to use this, fix later)
  if (strcmp(in.stls_guess_file,"NO_FILE")==0){
    for (int ii=0; ii < in.nx; ii++) {
      GG[ii] = 0.0; // Static local field correction
      GG_new[ii] = 1.0;
    }
    compute_ssf_stls(SS, SSHF, GG, phi, xx, in); // Static structure factor
  }
  else {
    read_guess_stls(SS, GG, in); // Read from file
  }

  // #######################################################################################
  // CALCULATIONS FOR A GIVEN VALUE OF ALPHA
  // #######################################################################################
  drs = rsArray[1] - rsArray[0];
  // Internal energy
  for (int ii=0; ii<nrs; ii++) {

    // State point
    in.rs = rsArray[ii];
   
    // Initialize arrays that depend only on the state point
    init_state_point_vs_stls_arrays(&in, xx, phi, SSHF, false);

    // Solve state point
    vs_stls_iterations(SS, SSHF, GG, GG_new, phi, xx, drs, alpha, in, false);
    
    // Internal energy
    uint[ii] = compute_internal_energy(SS, xx, in);
    printf("%f %f\n", in.rs, in.rs*uint[ii]);
    
  }

  // Free energy
  in.rs = rs;
  printf("Free energy: %f\n", compute_free_energy(uint, rsArray, in, nrs));


  // #######################################################################################
  // CALCULATIONS FOR A GIVEN VALUE OF ALPHA
  // #######################################################################################

  
  /* // Internal energy */
  /* if (verbose) printf("Internal energy: %.10f\n",compute_internal_energy(SS, xx, in)); */
  
  /* // Output to file */
  /* if (verbose) printf("Writing output files...\n"); */
  /* write_text_stls(SS, GG, phi, SSHF, xx, in); */
  /* write_guess_stls(SS, GG, in);  */
  /* if (verbose) printf("Done.\n"); */

  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
  free_vs_stls_arrays(uint, rsArray);
 
}


void alloc_vs_stls_arrays(input in, int nrs, double **uint, double **rsArray){
  *uint = malloc( sizeof(double) * nrs);
  *rsArray = malloc( sizeof(double) * nrs); 
}

void free_vs_stls_arrays(double *uint, double *rsArray){
  free(uint);
  free(rsArray);
}

void init_fixed_vs_stls_arrays(input *in, int nrs,  double *xx, double *rsArray, bool verbose){

  // Print on screen the parameter used to solve the STLS equation
  printf("------ Parameters used in the solution -------------\n");
  printf("Quantum degeneracy parameter: %f\n", in->Theta);
  printf("Quantum coupling parameter: %f\n", in->rs);
  printf("Chemical potential (low and high bound): %f %f\n", 
	 in->mu_lo, in->mu_hi);
  printf("Wave-vector cutoff: %f\n", in->xmax);
  printf("Wave-vector resolutions: %f\n", in->dx);
  printf("Number of Matsubara frequencies: %d\n", in->nl);
  printf("Maximum number of iterations: %d\n", in->nIter);
  printf("Error for convergence: %.5e\n", in->err_min_iter);
  printf("Number of threads: %d\n", omp_get_max_threads());
  printf("----------------------------------------------------\n");
  
  // Wave-vector grid
  if (verbose) printf("Wave-vector grid initialization: ");
  wave_vector_grid(xx, in);
  if (verbose) printf("Done.\n");

  // Array of coupling parameters
  if (verbose) printf("Coupling parameter grid initialization: ");
  rs_grid(rsArray, nrs, in);
  if (verbose) printf("Done.\n");
  
}


void init_state_point_vs_stls_arrays(input *in, double *xx,
				     double *phi, double *SSHF,
				     bool verbose){
 
  // Chemical potential
  if (in->Theta > 0) {
    if (verbose) printf("Chemical potential calculation: ");
    in->mu = compute_chemical_potential(*in);
    if (verbose) printf("Done. Chemical potential: %.8f\n", in->mu);
  }
  
  // Normalized ideal Lindhard density response
  if (in->Theta > 0) {
    if (verbose) printf("Normalized ideal Lindhard density calculation:\n");
    compute_idr(phi, xx, *in, verbose);
    if (verbose) printf("Done.\n");
  }
  
  // Static structure factor in the Hartree-Fock approximation
  if (verbose) printf("Static structure factor in the Hartree-Fock approximation: ");
  if (in->Theta == 0) {
    compute_ssf_HF_zero_temperature(SSHF, xx, *in);
  }
  else {
    compute_ssf_HF_finite_temperature(SSHF, xx, *in);
  }
  if (verbose) printf("Done.\n");

}


// -------------------------------------------------------------------
// FUNCTION USED TO PERFORM THE ITERATIONS FOR THE STLS SCHEME
// -------------------------------------------------------------------

void vs_stls_iterations(double *SS, double *SSHF,
			double *GG, double *GG_new,
			double *phi, double *xx,
			double drs, double alpha,
			input in, bool verbose) {

  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    double tic = omp_get_wtime();
    
    // Update SSF
    compute_ssf_stls(SS, SSHF, GG, phi, xx, in);

    // Update SLFC    
    compute_vs_slfc(GG_new, SS, xx, drs, alpha, in);

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
 
 
}


void compute_vs_slfc(double *GG, double *SS, double *xx, double drs, double alpha, input in) {

  // Arrays for VS-STLS solution
  double *GGrs = NULL;
  double *GGrsp = NULL;
  double *GGrsm = NULL;
  double rs = in.rs;

  // Allocate arrays
  GGrs = malloc( sizeof(double) * in.nx);
  GGrsp = malloc( sizeof(double) * in.nx);
  GGrsm = malloc( sizeof(double) * in.nx);
  
  // STLS local field correction at three state points
  compute_slfc(GGrs, SS, xx, in);
  in.rs = rs + drs;
  compute_slfc(GGrsp, SS, xx, in);
  in.rs = rs - drs;
  compute_slfc(GGrsm, SS, xx, in);

  // VS-STLS static local field correction
  for (int ii=0; ii<in.nx; ii++){

    // STLS contribution
    GG[ii] = GGrs[ii];

    // State point derivative contribution
    GG[ii] += -alpha/3.0*rs*(GGrsp[ii] - GGrsm[ii])/(2.0*drs);

    // Wave-vector derivative contribution
    if (ii > 0 && ii < in.nx-1) {
      GG[ii] += -alpha/3.0*xx[ii]*(GGrs[ii+1] - GGrs[ii-1])/(2.0*in.dx);
    }
    
  }
  
  // Free memory
  free(GGrs);
  free(GGrsp);
  free(GGrsm);
  
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FREE ENERGY
// -------------------------------------------------------------------

struct fex_params {

  gsl_spline *uint_sp_ptr;
  gsl_interp_accel *uint_acc_ptr;

};

double compute_free_energy(double *uint, double *rsArray, input in, int nrs) {

  double err;
  size_t neval;
  double fre;

  // Declare accelerator and spline objects
  gsl_spline *uint_sp_ptr;
  gsl_interp_accel *uint_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  uint_sp_ptr = gsl_spline_alloc(gsl_interp_linear, nrs);
  uint_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(uint_sp_ptr, rsArray, uint, nrs);

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &fex;

  // Internal energy
  struct fex_params fexp = {uint_sp_ptr, uint_acc_ptr};
  ff_int.params = &fexp;
  gsl_integration_cquad(&ff_int,
  			rsArray[0], in.rs,
  			0.0, QUAD_REL_ERR,
  			wsp,
  			&fre, &err, &neval);
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(uint_sp_ptr);
  gsl_interp_accel_free(uint_acc_ptr);

  // Output
  return fre/(in.rs*in.rs);

}

double fex(double rs, void* pp) {

  struct fex_params* params = (struct fex_params*)pp;
  gsl_spline* uint_sp_ptr = (params->uint_sp_ptr);
  gsl_interp_accel* uint_acc_ptr = (params->uint_acc_ptr);

  return rs*gsl_spline_eval(uint_sp_ptr, rs, uint_acc_ptr);
  
}


// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rsArray, int nrs, input *in){

  double  drs = in->rs/nrs;
  for (int ii=0; ii<nrs; ii++) {
    rsArray[ii] = drs + drs*(ii);
  }
  
}
