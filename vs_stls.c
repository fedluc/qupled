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
// FUNCTION USED TO SOLVE THE VS-STLS SCHEME
// -------------------------------------------------------------------

void solve_vs_stls(input in, bool verbose) {

  // Arrays for VS-STLS solution
  vs_sp xx;
  vs_sp phi;
  vs_sp GG;
  vs_sp GG_new;
  vs_sp SS;
  vs_sp SSHF;
  vs_sp rsu;
  vs_sp rsa;
  input vs_in[VS_SP_EL];

  // Limit to zero temperature
  /* if (in.Theta > 0.0) { */
  /*   printf("The VS-STLS scheme is only implemented in the ground state\n"); */
  /*   exit(EXIT_FAILURE); */
  /* } */
  
  // Allocate arrays
  stls_pointers stls_pp;
  vs_stls_pointers vs_stls_pp;
  for (int ii=0; ii<VS_SP_EL; ii++) {
    alloc_stls_arrays_new(in, &stls_pp);
    alloc_vs_stls_struct_arrays(in, stls_pp, &xx,
    				&phi, &SS, &SSHF,
    				&GG, &GG_new, ii);
    alloc_vs_stls_thermo_arrays_elem(in, &vs_stls_pp);
    alloc_vs_stls_thermo_arrays(in, vs_stls_pp, &rsu, &rsa, ii);
  }
  
  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_vs_stls_arrays(&in, vs_in, xx, rsa, verbose);
  init_state_point_vs_stls_arrays(vs_in, xx, phi, SSHF, verbose);
   
  // Iterative procedure to fix the compressibility sum rule
  if (verbose) printf("Iterative procedure to enforce the compressibility sum rule ...\n");
  in.a_csr = vs_stls_thermo_iterations(xx, rsu, rsa, vs_in, true);
  if (verbose) printf("Done.\n");

  // Structural properties
  if (verbose) printf("Structural properties calculation...\n");
  initial_guess_vs_stls(xx, SS, SSHF, GG, GG_new, phi, vs_in);
  vs_stls_struct_iterations(SS, SSHF, GG, GG_new, phi, xx, vs_in, false);
  if (verbose) printf("Done.\n");
  
  // Thermodynamic properties
  if (verbose) printf("Free parameter: %.5f\n",vs_in[0].a_csr);
  if (verbose) printf("Internal energy: %.10f\n",compute_internal_energy(SS.in, xx.in, in));
  if (verbose) printf("Free energy: %.10f\n",compute_free_energy(rsu.in, rsa.in, in));
  
  // Output to file
  /* if (verbose) printf("Writing output files...\n"); */
  /* write_text_stls(SS.in, GG.in, phi.in, SSHF.in, xx.in, in); */
  /* write_text_vs_stls(rsu, rsa, in); */
  /* write_guess_stls(SS.in, GG.in, in); */
  /* if (verbose) printf("Done.\n"); */

  // Free memory
  /* free_vs_stls_struct_arrays(xx, phi, SS, SSHF, GG, GG_new); */
  /* free_vs_stls_thermo_arrays(rsu, rsa); */
  
}


// -------------------------------------------------------------------
// FUNCTION USED TO LOOP OVER THE VS_SP DATA STRUCTURES
// -------------------------------------------------------------------

double *get_el(vs_sp sp, int ii) {
  
  switch(ii) {
    
  case 0:
    return sp.in;
  case 1:
    return sp.rsp1;
  case 2:
    return sp.rsm1;
  case 3:
    return sp.rsp2;
  case 4:
    return sp.rsm2;
  case 5:
    return sp.tp1;
  case 6:
    return sp.tm1;
  case 7:
    return sp.tp2;
  case 8:
    return sp.tm2;
    
  }

  printf("Invalid structure element for vs_sp data structures\n"); 
  exit(EXIT_FAILURE);

}

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

// Function to allocate the VS-STLS arrays for the structural properties
void alloc_vs_stls_struct_arrays(input in, stls_pointers pp, vs_sp *xx,
				 vs_sp *phi, vs_sp *SS, vs_sp *SSHF,
				 vs_sp *GG, vs_sp *GG_new, int el){

  bool alloc_failed = false;
  
  if (pp.xx == NULL || pp.phi == NULL || pp.SS == NULL ||
      pp.SSHF == NULL || pp.GG == NULL || pp.GG_new == NULL)
    alloc_failed = true;
  
  switch(el) {
    
  case 0:
    xx->in = pp.xx;    
    phi->in = pp.phi;
    SS->in = pp.SS;
    SSHF->in = pp.SSHF;
    GG->in = pp.GG;
    GG_new->in= pp.GG_new;
    break;
    
  case 1:
    xx->rsp1 = pp.xx;    
    phi->rsp1 = pp.phi;
    SS->rsp1 = pp.SS;
    SSHF->rsp1 = pp.SSHF;
    GG->rsp1 = pp.GG;
    GG_new->rsp1 = pp.GG_new;
    break;

  case 2:
    xx->rsm1 = pp.xx;    
    phi->rsm1 = pp.phi;
    SS->rsm1 = pp.SS;
    SSHF->rsm1 = pp.SSHF;
    GG->rsm1 = pp.GG;
    GG_new->rsm1 = pp.GG_new;
    break;

  case 3:
    xx->rsp2 = pp.xx;    
    phi->rsp2 = pp.phi;
    SS->rsp2 = pp.SS;
    SSHF->rsp2 = pp.SSHF;
    GG->rsp2 = pp.GG;
    GG_new->rsp2 = pp.GG_new;
    break;

  case 4:
    xx->rsm2 = pp.xx;    
    phi->rsm2 = pp.phi;
    SS->rsm2 = pp.SS;
    SSHF->rsm2 = pp.SSHF;
    GG->rsm2 = pp.GG;
    GG_new->rsm2 = pp.GG_new;
    break;

  case 5:
    xx->tp1 = pp.xx;    
    phi->tp1 = pp.phi;
    SS->tp1 = pp.SS;
    SSHF->tp1 = pp.SSHF;
    GG->tp1 = pp.GG;
    GG_new->tp1 = pp.GG_new;
    break;

  case 6:
    xx->tp2 = pp.xx;    
    phi->tp2 = pp.phi;
    SS->tp2 = pp.SS;
    SSHF->tp2 = pp.SSHF;
    GG->tp2 = pp.GG;
    GG_new->tp2 = pp.GG_new;
    break;

  case 7:
    xx->tm1 = pp.xx;    
    phi->tm1 = pp.phi;
    SS->tm1 = pp.SS;
    SSHF->tm1 = pp.SSHF;
    GG->tm1 = pp.GG;
    GG_new->tm1 = pp.GG_new;
    break;

  case 8:
    xx->tm2 = pp.xx;    
    phi->tm2 = pp.phi;
    SS->tm2 = pp.SS;
    SSHF->tm2 = pp.SSHF;
    GG->tm2 = pp.GG;
    GG_new->tm2 = pp.GG_new;
    break;

  default:
    alloc_failed = true;

  }

  if (alloc_failed) {
    printf("The allocation of the VS-STLS structural arrays failed\n"); 
    exit(EXIT_FAILURE);
  }
 
}

// Function to allocate the VS-STLS arrays for the thermodynamic properties
void alloc_vs_stls_thermo_arrays(input in, vs_stls_pointers pp,
				 vs_sp *rsu, vs_sp *rsa, int el){

  bool alloc_failed = false;
    
  if (pp.rsu == NULL ||
      pp.rsa == NULL)
    alloc_failed = true;
  
  switch(el) {
    
  case 0:
    rsu->in = pp.rsu;
    rsa->in = pp.rsa;
    break;
    
  case 1:
    rsu->rsp1 = pp.rsu;
    rsa->rsp1 = pp.rsa;
    break;

  case 2:
    rsu->rsm1 = pp.rsu;
    rsa->rsm1 = pp.rsa;
    break;

  case 3:
    rsu->rsp2 = pp.rsu;
    rsa->rsp2 = pp.rsa;
    break;

  case 4:
    rsu->rsm2 = pp.rsu;
    rsa->rsm2 = pp.rsa;
    break;

  case 5:
    rsu->tp1 = pp.rsu;
    rsa->tp1 = pp.rsa;
    break;

  case 6:
    rsu->tm1 = pp.rsu;
    rsa->tm1 = pp.rsa;
    break;

  case 7:
    rsu->tp2 = pp.rsu;
    rsa->tp2 = pp.rsa;
    break;

  case 8:
    rsu->tm2 = pp.rsu;
    rsa->tm2 = pp.rsa;
    break;

    
  default:
    alloc_failed = true;

  }

  if (alloc_failed) {
    printf("The allocation of the VS-STLS thermodynamic arrays failed\n"); 
    exit(EXIT_FAILURE);
  }
  
}

// Function to allocate one element of the VS-STLS arrays for the thermodynamic properties
void alloc_vs_stls_thermo_arrays_elem(input in, vs_stls_pointers *pp){

  pp->rsu = malloc( sizeof(double) * in.nrs);
  if (pp->rsu == NULL) {
    fprintf(stderr, "Failed to allocate memory for the free energy integrand\n");
    exit(EXIT_FAILURE);
  }

  pp->rsa = malloc( sizeof(double) * in.nrs);
  if (pp->rsa == NULL) {
    fprintf(stderr, "Failed to allocate memory for the coupling parameter array\n");
    exit(EXIT_FAILURE);
  }

  
}


// Function to allocate a single VS-STLS array
void alloc_vs_stls_one_array(double *pp, vs_sp *vs_arr, int el){

  bool alloc_failed = false;

  // Allocate vs_sp array
  if (pp == NULL)
    alloc_failed = true;
  
  switch(el) {
    
  case 0:
    vs_arr->in = pp;
    break;
    
  case 1:
    vs_arr->rsp1 = pp;
    break;

  case 2:
    vs_arr->rsm1 = pp;
    break;

  case 3:
    vs_arr->rsp2 = pp;
    break;

  case 4:
    vs_arr->rsm2 = pp;
    break;

  case 5:
    vs_arr->tp1 = pp;
    break;

  case 6:
    vs_arr->tm1 = pp;
    break;

  case 7:
    vs_arr->tp2 = pp;
    break;

  case 8:
    vs_arr->tm2 = pp;
    break;

  default:
    alloc_failed = true;

  }

  if (alloc_failed) {
    printf("The allocation of the VS-STLS array failed\n"); 
    exit(EXIT_FAILURE);
  }
 
}

// Function to free the VS-STLS arrays for the structural properties
void free_vs_stls_struct_arrays(vs_sp xx, vs_sp phi, vs_sp SS, vs_sp SSHF,
				vs_sp GG, vs_sp GG_new){

  for (int ii=0; ii<VS_SP_EL; ii++) {
    free_stls_arrays(get_el(xx,ii), get_el(phi,ii), get_el(GG,ii),
		     get_el(GG_new, ii), get_el(SS,ii),
		     get_el(SSHF,ii));
  }
 
}

// Function to free the VS-STLS arrays for the thermodynamic properties
void free_vs_stls_thermo_arrays(vs_sp rsu, vs_sp rsa){

  for (int ii=0; ii<VS_SP_EL; ii++) {
    free(get_el(rsu,ii));
    free(get_el(rsa,ii));
  }
 
}


// Function to free a single VS-STLS array
void free_vs_stls_one_array(vs_sp vs_arr){

  for (int ii=0; ii<VS_SP_EL; ii++) {
    free(get_el(vs_arr,ii));
  }
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

// Initialize arrays that do not depend on iterations and state points
void init_fixed_vs_stls_arrays(input *in, input *vs_in,
			       vs_sp xx, vs_sp rsa,
			       bool verbose){

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
  for (int ii=0; ii<VS_SP_EL; ii++) {
    wave_vector_grid(get_el(xx,ii), in);
  }
  if (verbose) printf("Done.\n");

  // Array of coupling parameters
  if (verbose) printf("Coupling parameter grid initialization: ");
  for (int ii=0; ii<VS_SP_EL; ii++) {
    rs_grid(get_el(rsa,ii), in);
  }
  if (verbose) printf("Done.\n");

  // Input structure for the state point of interest
  vs_in[0] = *in;
  
}


// Initialize arrays that do not depend on the iterations, but that
// are a function of the state point
void init_state_point_vs_stls_arrays(input *vs_in, vs_sp xx,
				     vs_sp phi, vs_sp SSHF,
				     bool verbose){

  input in = vs_in[0];

  // Input structure for points to be solved simultaneously
  for(int ii=1; ii<VS_SP_EL; ii++){
    vs_in[ii] = in;			       
  }
  vs_in[1].rs += in.drs;
  vs_in[2].rs -= in.drs;
  vs_in[3].rs += 2.0*in.drs;
  vs_in[4].rs -= 2.0*in.drs;
  if (vs_in[2].rs < 0.0) vs_in[2].rs = 0.0;
  if (vs_in[4].rs < 0.0) vs_in[4].rs = 0.0;
  if (VS_SP_EL > 5) {
    vs_in[5].Theta += in.drs;
    vs_in[6].Theta -= in.drs;
    vs_in[7].Theta += 2.0*in.drs;
    vs_in[8].Theta -= 2.0*in.drs;
    if (vs_in[7].rs < 0.0) vs_in[7].rs = 0.0;
    if (vs_in[8].rs < 0.0) vs_in[8].rs = 0.0;
  }
  
  // Chemical potential
  if (verbose) printf("Chemical potential calculation: ");
  compute_chemical_potential_vs_stls(vs_in);
  if (verbose) printf("Done. Chemical potential: %.8f\n", in.mu);
  
  // Normalized ideal Lindhard density response
  if (verbose) printf("Normalized ideal Lindhard density calculation: ");
  compute_idr_vs_stls(phi, xx, vs_in);
  if (verbose) printf("Done.\n");
  
  // Static structure factor in the Hartree-Fock approximation
  if (verbose) printf("Static structure factor in the Hartree-Fock approximation: ");
  compute_ssf_HF_vs_stls(SSHF, xx, vs_in);
  if (verbose) printf("Done.\n");

}


// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rsa, input *in){

  rsa[0] = 0.0;
  for (int ii=1; ii<in->nrs; ii++) {
    rsa[ii] =  in->drs*(ii);
  }
  
}

// ---------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE INITIAL GUESS
// ---------------------------------------------------------------------
void initial_guess_vs_stls(vs_sp xx, vs_sp SS, vs_sp SSHF,
			   vs_sp GG, vs_sp GG_new, vs_sp phi,
			   input *vs_in){

  input in = vs_in[0];
  
  if (strcmp(in.stls_guess_file,"NO_FILE")==0){

    for (int ii=0; ii<VS_SP_EL; ii++){
      initial_guess_stls(get_el(xx,ii), get_el(SS,ii), get_el(SSHF, ii),
			 get_el(GG,ii), get_el(GG_new,ii), get_el(phi,ii),
			 vs_in[ii]);
    }

  }  
  else {

    // Read from file
    read_guess_vs_stls(SS, GG, vs_in);
    
  }
  
}

// ---------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE VS-STLS SCHEME
// ---------------------------------------------------------------------

// Iterations over the parameter used to enforce the CSR rule
double vs_stls_thermo_iterations(vs_sp xx, vs_sp rsu, vs_sp rsa,
				 input *vs_in, bool verbose) {

  input in = vs_in[0];
  double alpha = in.a_csr;

  // Iterations
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {

    // Start timing
    double tic = omp_get_wtime();

    // Get parameter to enforce the compressibility sum rule
    alpha = compute_alpha(xx, rsu, rsa, vs_in);

    // Update diagnostic
    iter_counter++;
    iter_err = vs_stls_thermo_err(alpha, vs_in);
    vs_stls_thermo_update(alpha, vs_in);

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

  
  return alpha;
  
}

double vs_stls_thermo_err(double alpha, input *vs_in){

  return fabs((alpha - vs_in[0].a_csr)/alpha);
  
}

void vs_stls_thermo_update(double alpha, input *vs_in){

  for (int ii=0; ii<VS_SP_EL; ii++) {
    vs_in[ii].a_csr = alpha;
  }
   
}


// Iterations over the structural properties
void vs_stls_struct_iterations(vs_sp SS, vs_sp SSHF,
			       vs_sp GG, vs_sp GG_new,
			       vs_sp phi, vs_sp xx,
			       input *vs_in, bool verbose) {

  input in = vs_in[0];
    
  // Iterations
  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter) {
    
    // Start timing
    double tic = omp_get_wtime();
    
    // Update SSF
    compute_ssf_vs_stls(SS, SSHF, GG, phi, xx, vs_in);
    
    // Update SLFC
    compute_slfc_vs_stls(GG_new, SS, xx, vs_in);
    
    // Update diagnostic
    iter_counter++;
    iter_err = vs_stls_struct_err(GG, GG_new, vs_in);
    vs_stls_struct_update(GG, GG_new, vs_in);
    
    // End timing
    double toc = omp_get_wtime();
    
    // Print diagnostic
    if (verbose) {
      printf("--- iteration (structure) %d ---\n", iter_counter);
      printf("Elapsed time: %f seconds\n", toc - tic);
      printf("Residual error: %.5e\n", iter_err);
      fflush(stdout);
    }
    
  }

  // Stop execution if convergence was not reached
  /* if (iter_err > in.err_min_iter) { */
  /*   fprintf(stderr, "The calculation for the structural properties for " */
  /* 	    "the state point (rs = %f, theta = %f) did not converge\n", */
  /* 	    vs_in[0].rs, vs_in[0].Theta); */
  /*   exit(EXIT_FAILURE); */
  /* } */
	    
  if (verbose) printf("Done.\n");
 
}

double vs_stls_struct_err(vs_sp GG, vs_sp GG_new, input *vs_in){

  double err = 0.0;
  double err_tmp;
  
  for (int ii=0; ii<VS_SP_EL; ii++) {
    err_tmp = stls_err(get_el(GG,ii), get_el(GG_new,ii), vs_in[ii]);
    if (err_tmp > err) err = err_tmp;
  }
  return err;
  
}

void vs_stls_struct_update(vs_sp GG, vs_sp GG_new, input *vs_in){

  for (int ii=0; ii<VS_SP_EL; ii++) {
    stls_update(get_el(GG,ii), get_el(GG_new,ii), vs_in[ii]);
  }
   
}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

void compute_chemical_potential_vs_stls(input *vs_in){

  for (int ii=0; ii<VS_SP_EL; ii++){
    if (vs_in[ii].Theta > 0.0) {
      vs_in[ii].mu = compute_chemical_potential(vs_in[ii]);
    }
  }
  
}
					
// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_vs_stls(vs_sp SS, vs_sp SSHF, vs_sp GG, vs_sp phi,
		    vs_sp xx, input *vs_in){

  for (int ii=0; ii<VS_SP_EL; ii++) {
    compute_ssf_stls(get_el(SS,ii), get_el(SSHF,ii),
		     get_el(GG,ii), get_el(phi,ii),
		     get_el(xx,ii), vs_in[ii]);
  }
  
}

void compute_ssf_HF_vs_stls(vs_sp SS, vs_sp xx, input *vs_in){

  for (int ii=0; ii<VS_SP_EL; ii++) {
    if (vs_in[ii].Theta == 0) {
      compute_ssf_HF_zero_temperature(get_el(SS,ii), get_el(xx,ii),
				      vs_in[ii]);
    }
    else {
      compute_ssf_HF_finite_temperature(get_el(SS,ii), get_el(xx,ii),
					vs_in[ii]);
    }
	
  }
  
}
		    
// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_vs_stls(vs_sp GG, vs_sp SS, vs_sp xx, input *vs_in) {

  bool finite_temperature = false;
  vs_sp GG_stls;
  double *pp;
  input in = vs_in[0];
  double alpha = in.a_csr;
  double a_drs = alpha/(3.0*2.0*in.drs);
  double a_dx = alpha/(3.0*2.0*in.dx);
  double a_dt = 2.0*alpha/(3.0*2.0*in.drs);
  /* double der_coeff = 2.0; */
  
  // Check if finite temperature calculations must be performed
  if (in.Theta > 0.0) finite_temperature = true;
  
  // Allocate arrays
  for (int ii=0; ii<VS_SP_EL; ii++){
    pp = malloc(sizeof(double) * in.nx);
    alloc_vs_stls_one_array(pp, &GG_stls, ii);
  }

  
  // STLS contribution at all the state points that have to be updated simultaneously
  for (int ii=0; ii<VS_SP_EL; ii++) {
    compute_slfc(get_el(GG_stls,ii), get_el(SS,ii), get_el(xx, ii), vs_in[ii]);
  }

  // STLS contribution to the VS-STLS static local field correction
  for (int ii=0; ii<VS_SP_EL; ii++){
    for (int jj=0; jj<in.nx; jj++){
      get_el(GG, ii)[jj] = get_el(GG_stls,ii)[jj];
    }
  }

  // Wave-vector derivative contribution to the VS-STLS static local field correction
  for (int ii=0; ii<VS_SP_EL; ii++){
    for (int jj=1; jj<in.nx-1; jj++){
      get_el(GG,ii)[jj] += -a_dx*get_el(xx,ii)[jj]*(get_el(GG_stls,ii)[jj+1]
						    - get_el(GG_stls,ii)[jj-1]);
    }
  }

  // State point derivative contribution to the VS-STLS static local field correction
  for (int ii=0; ii<in.nx; ii++) {
    
    // State point derivative contribution
    GG.in[ii] += -a_drs*vs_in[0].rs
                *(GG_stls.rsp1[ii] - GG_stls.rsm1[ii]);
    GG.rsp1[ii] += -a_drs*vs_in[1].rs
                  *(GG_stls.rsp2[ii] - GG_stls.in[ii]);
    GG.rsm1[ii] += -a_drs*vs_in[2].rs
                   *(GG_stls.in[ii] - GG_stls.rsm2[ii]);
    GG.rsp2[ii] += -a_drs*vs_in[3].rs
                  *(3.0*GG_stls.rsp2[ii] - 4.0*GG_stls.rsp1[ii] + GG_stls.in[ii]);
    GG.rsm2[ii] += -a_drs*vs_in[4].rs
                  *(-3.0*GG_stls.rsm2[ii] + 4.0*GG_stls.rsm1[ii] - GG_stls.in[ii]);
    if (finite_temperature) {
      GG.in[ii] += -a_dt*vs_in[0].Theta
    	           *(GG_stls.tp1[ii] - GG_stls.tm1[ii]);
      GG.tp1[ii] += -a_dt*vs_in[5].Theta
                   *(GG_stls.tp2[ii] - GG_stls.in[ii]);
      GG.tm1[ii] += -a_dt*vs_in[6].Theta
                    *(GG_stls.in[ii] - GG_stls.tm2[ii]);
      GG.tp2[ii] += -a_dt*vs_in[7].Theta
    	            *(3.0*GG_stls.tp2[ii] - 4.0*GG_stls.tp1[ii] + GG_stls.in[ii]);
      GG.tm2[ii] += -a_dt*vs_in[8].Theta
                     *(-3.0*GG_stls.tm2[ii] + 4.0*GG_stls.tm1[ii] - GG_stls.in[ii]);
      
    }

  }
  
  // Free memory
  free_vs_stls_one_array(GG_stls);
    
}

// -----------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY RESPONSE
// -----------------------------------------------------------------------

void compute_idr_vs_stls(vs_sp phi, vs_sp xx, input *vs_in){

  for (int ii=0; ii<VS_SP_EL; ii++) {
    if (vs_in[ii].Theta > 0.0) {
      compute_idr(get_el(phi,ii), get_el(xx,ii),
		  vs_in[ii], false);
    }
  }
  
}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE PARAMETER FOR THE COMPRESSIBILITY RULE
// -------------------------------------------------------------------

double compute_alpha(vs_sp xx, vs_sp rsu, vs_sp rsa, input *vs_in){

  bool finite_temperature = false;
  double urs = 0.0, ursp = 0.0, ursm = 0.0;
  double utp = 0.0, utm = 0.0;
  double frs = 0.0, frsp = 0.0, frsm = 0.0;
  double ftp = 0.0, ftm = 0.0;
  double frsptp = 0.0, frsmtp = 0.0;
  double frsptm = 0.0, frsmtm = 0.0;
  double dudrs, dudt;
  double dfdrs, dfdt;
  double d2fdrs2, d2fdt2, d2fdrsdt;
  double numer;
  double denom;
  double alpha;
  input in = vs_in[0];
  int nrs = in.nrs;
  
  // Check if finite temperature calculations must be performed
  if (in.Theta > 0.0) finite_temperature = true;
  
  // Get integrand for the free energy
  compute_rsu(xx, rsu, rsa, vs_in, false);

  // Internal energy
  urs = rsu.in[nrs-2]/rsa.in[nrs-2];
  ursp = rsu.in[nrs-1]/rsa.in[nrs-1];
  ursm = rsu.in[nrs-3]/rsa.in[nrs-3];
  if (finite_temperature) {
    utp = rsu.tp1[nrs-2]/rsa.tp1[nrs-2];
    utm = rsu.tm1[nrs-2]/rsa.tm1[nrs-2];
  } 
  
  // Free energy
  frs = compute_free_energy(rsu.in, rsa.in, vs_in[0]);
  frsp = compute_free_energy(rsu.in, rsa.in, vs_in[1]);
  frsm = compute_free_energy(rsu.in, rsa.in, vs_in[2]);
  if (finite_temperature) {
    ftp = compute_free_energy(rsu.tp1, rsa.tp1, vs_in[0]);
    ftm = compute_free_energy(rsu.tm1, rsa.tm1, vs_in[0]);
    frsptp = compute_free_energy(rsu.tp1, rsa.tp1, vs_in[1]);
    frsptm = compute_free_energy(rsu.tm1, rsa.tm1, vs_in[1]);
    frsmtp = compute_free_energy(rsu.tp1, rsa.tp1, vs_in[2]);
    frsmtp = compute_free_energy(rsu.tm1, rsa.tm1, vs_in[2]);
  }
  
  // Internal energy derivatives
  dudrs = (ursp - ursm)/(2.0*in.drs);
  if (finite_temperature)
    dudt =  (utp - utm)/(2.0*in.drs);
  
  // Free energy derivatives
  dfdrs = (frsp - frsm)/(2.0*in.drs);
  d2fdrs2 = (frsp -2.0*frs + frsm)/(in.drs*in.drs);
  if (finite_temperature) {
    dfdt = (ftp - ftm)/(2.0*in.drs);
    d2fdt2 = (ftp -2.0*frs + ftm)/(in.drs*in.drs);
    d2fdrsdt = (frsptp - frsmtp - frsptm + frsmtm)/(4.0*in.drs*in.drs);
  }
    
  // Parameter for the compressibility sum rule
  numer = (2.0*frs - (1.0/6.0)*in.rs*in.rs*d2fdrs2
  	   + (4.0/3.0)*in.rs*dfdrs);
  denom = (urs + (1.0/3.0)*in.rs*dudrs);
  if (finite_temperature) {
    numer += -(2.0/3.0)*in.Theta*in.Theta*d2fdt2
             -(2.0/3.0)*in.Theta*in.rs*d2fdrsdt
             +(1.0/3.0)*in.Theta*dfdt;
    denom += (2.0/3.0)*in.Theta*dudt;
  }
  alpha = numer/denom;

  // Output
  return alpha;
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE INTEGRAND FOR THE FREE ENERGY
// -------------------------------------------------------------------

void compute_rsu(vs_sp xx, vs_sp rsu, vs_sp rsa,
		 input *vs_in, bool verbose){

  int nrs = vs_in[0].nrs;
  
  #pragma omp parallel
  {
    #pragma omp for // Parallel calculations
    for (int ii=0; ii<nrs; ii++) {
 
      // Arrays for VS-STLS solution
      vs_sp tmp_xx;
      vs_sp phi;
      vs_sp GG;
      vs_sp GG_new;
      vs_sp SS;
      vs_sp SSHF;
      input vs_in_tmp[VS_SP_EL];
      double u_int;
      stls_pointers stls_pp;

      // Define state point
      for (int jj=0; jj<VS_SP_EL; jj++) {
	vs_in_tmp[jj] = vs_in[jj];
      }
      vs_in_tmp[0].rs = rsa.in[ii];
      
      // Allocate arrays
      for (int jj=0; jj<VS_SP_EL; jj++) {
	alloc_stls_arrays_new(vs_in[0], &stls_pp);
	alloc_vs_stls_struct_arrays(vs_in[0], stls_pp, &tmp_xx, &phi,
				    &SS, &SSHF, &GG, &GG_new, jj);
      }

      // Initialize arrays that depend only on the state point
      init_state_point_vs_stls_arrays(vs_in_tmp, xx, phi, SSHF, verbose);
  
      // Initial guess (same for all threads)
      // NOTE: This could be problematic if the initial
      // rs is noticeably different from zero (rs > 20)
      initial_guess_vs_stls(xx, SS, SSHF, GG, GG_new, phi, vs_in_tmp);
      
      // Compute structural properties with the VS-STLS static local field correction
      vs_stls_struct_iterations(SS, SSHF, GG, GG_new, phi, xx, vs_in_tmp, verbose);
      
      // Compute the free energy integrand
      for (int jj=0; jj<VS_SP_EL; jj++){

	// Avoid divergencies for rs = 0.0 (the internal energy scales as 1.0/rs)
	if (vs_in_tmp[jj].rs == 0.0) vs_in_tmp[jj].rs = 1.0; // Avoid divergencies for rs = 0.0;
	u_int = compute_internal_energy(get_el(SS,jj), get_el(xx,jj), vs_in_tmp[jj]);
	get_el(rsu,jj)[ii] = vs_in_tmp[jj].rs*u_int;
	
      }
      
      // Free memory
      free_vs_stls_struct_arrays(tmp_xx, phi, SS, SSHF, GG, GG_new);
      
    }

 }
  
}



// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

// write text files for output
void write_text_vs_stls(double *rsu, double *rsa, input in){

    FILE* fid;
    char out_name[100];
    
    // Output for the internal energy used for thermodynamic integration
    sprintf(out_name, "rsu_thermoint_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        fprintf(stderr, "Error while creating the output file for the static structure factor (HF)");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nrs; ii++)
        fprintf(fid, "%.8e %.8e\n", rsa[ii], rsu[ii]);

    fclose(fid);

    // Output for the free energy
    sprintf(out_name, "fxc_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        fprintf(stderr, "Error while creating the output file for the free energy\n");
        exit(EXIT_FAILURE);
    }
    fprintf(fid, "%.8e\n", compute_free_energy(rsu, rsa, in));
    fclose(fid);

    // Output for the free parameter
    sprintf(out_name, "alpha_csr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        fprintf(stderr, "Error while creating the output file for the free parameter\n");
        exit(EXIT_FAILURE);
    }
    fprintf(fid, "%.8e\n", in.a_csr);
    fclose(fid);

}

// Read guess from input file
void read_guess_vs_stls(vs_sp SS, vs_sp GG, input *vs_in){
  
  for (int ii=0;  ii<VS_SP_EL; ii++) {
    read_guess_stls(get_el(SS,ii), get_el(GG,ii),
		    vs_in[ii]);
  }
  
}
