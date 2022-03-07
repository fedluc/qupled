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
  vs_struct xx;
  vs_struct phi;
  vs_struct GG;
  vs_struct GG_new;
  vs_struct SS;
  vs_struct SSHF;
  vs_struct rsu;
  vs_struct rsa;
  input vs_in[VSS_NUMEL];
  
  // Temporary arrays for allocation
  double *xx_el = NULL; 
  double *phi_el = NULL;
  double *GG_el = NULL;
  double *GG_new_el = NULL;
  double *SS_el = NULL;
  double *SSHF_el = NULL;
  double *rsu_el = NULL;
  double *rsa_el = NULL;

  // Define size of grid for thermodynamic integration
  get_grid_thermo_size(&in);
  
  // Allocate arrays
  for (int ii=0; ii<VSS_NUMEL; ii++) {
    alloc_stls_arrays(in, &xx_el, &phi_el, &GG_el,
		      &GG_new_el, &SS_el, &SSHF_el);
    alloc_vs_stls_arrays(in, &rsa_el, &rsu_el);
    xx.el[ii] = xx_el;
    phi.el[ii] = phi_el;
    GG.el[ii] = GG_el;
    GG_new.el[ii] = GG_new_el;
    SS.el[ii] = SS_el;
    SSHF.el[ii] = SSHF_el;
    rsa.el[ii] = rsa_el;
    rsu.el[ii] = rsu_el;
  }

  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_vs_stls_arrays(&in, vs_in, xx, rsa, verbose);
  init_tmp_vs_stls_arrays(vs_in, xx, phi, SSHF, verbose);

  // Iterative procedure to enforce the compressibility sum rule
  if (verbose) printf("Iterative procedure to enforce the compressibility sum rule ...\n");
  if (in.vs_solve_csr) 
    in.vs_alpha = vs_stls_thermo_iterations(xx, rsu, rsa, vs_in, true);
  if (verbose) printf("Done.\n");

  // Structural properties
  if (verbose) printf("Structural properties calculation...\n");
  initial_guess_vs_stls(xx, SS, SSHF, GG, GG_new, phi, vs_in);
  vs_stls_struct_iterations(SS, SSHF, GG, GG_new, phi, xx, vs_in, false);
  if (verbose) printf("Done.\n");
  
  // Thermodynamic properties
  if (verbose) printf("Free parameter: %.5f\n",vs_in[VSS_IDXIN].vs_alpha);
  if (verbose) printf("Internal energy: %.10f\n",compute_internal_energy(SS.rst, xx.rst, in));
  if (verbose) printf("Free energy: %.10f\n",compute_free_energy(rsu.rst, rsa.rst, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text_vs_stls(SS.rst, GG.rst, phi.rst, SSHF.rst, xx.rst,
  		     rsu.rst, rsa.rst, in);
  write_guess_stls(SS.rst, GG.rst, in);
  write_thermo_vs_stls(rsa, rsu, vs_in);
  if (verbose) printf("Done.\n");

  // Free memory
  for (int ii=0; ii<VSS_NUMEL; ii++){
    free_stls_arrays(xx.el[ii], phi.el[ii], GG.el[ii],
  		     GG_new.el[ii], SS.el[ii], SSHF.el[ii]);
    free_vs_stls_arrays(rsa.el[ii], rsu.el[ii]);
  }
  
}

// -------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE SIZE OF THE GRID FOR THERMODYNAMIC
// INTEGRATION
// -------------------------------------------------------------------

void get_grid_thermo_size(input *in){
  
  // Define a grid from 0.0 up to r_s + 2.0*d_rs  
  in->vs_nrs = (int)floor((in->rs + in->vs_drs)/in->vs_drs) + 2;
  
}

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

// Function to allocate one element of the VS-STLS arrays for the thermodynamic properties
void alloc_vs_stls_arrays(input in, double **rsa, double **rsu){

  *rsa = malloc( sizeof(double) * in.vs_nrs);
  if (*rsa == NULL) {
    fprintf(stderr, "Failed to allocate memory for the coupling parameter array\n");
    exit(EXIT_FAILURE);
  }

  *rsu = malloc( sizeof(double) * in.vs_nrs);
  if (*rsu == NULL) {
    fprintf(stderr, "Failed to allocate memory for the free energy integrand\n");
    exit(EXIT_FAILURE);
  }
  
}

// Function to free the VS-STLS arrays for the thermodynamic properties
void free_vs_stls_arrays(double *rsa, double *rsu){

  free(rsa);
  free(rsu);

}

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

// Initialize arrays that do not depend on iterations and state points
void init_fixed_vs_stls_arrays(input *in, input *vs_in,
			       vs_struct xx, vs_struct rsa,
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
  for (int ii=0; ii<VSS_NUMEL; ii++) {
    wave_vector_grid(xx.el[ii], in);
  }
  if (verbose) printf("Done.\n");

  // Array of coupling parameters
  if (verbose) printf("Coupling parameter grid initialization: ");
  for (int ii=0; ii<VSS_NUMEL; ii++) {
    rs_grid(rsa.el[ii], in);
  }
  if (verbose) printf("Done.\n");
  
  // Input structure for the state point of interest
  vs_in[VSS_IDXIN] = *in;
  
}


// Initialize arrays that do not depend on the iterations, but that
// are a function of the state point
void init_tmp_vs_stls_arrays(input *vs_in, vs_struct xx,
				     vs_struct phi, vs_struct SSHF,
				     bool verbose){

  bool finite_temperature = false;
  input in = vs_in[VSS_IDXIN];

  // Check if finite temperature calculations must be performed
  if (in.Theta > 0.0) finite_temperature = true;
  
  // Input structure for points to be solved simultaneously
  for(int ii=0; ii<VSS_NUMEL; ii++){
    vs_in[ii] = in;
  }
  
  // State points with modified coupling parameter
  for (int ii=0; ii<VSS_STENCIL; ii++){
    for(int jj=ii; jj<VSS_NUMEL; jj=jj+VSS_STENCIL){
      vs_in[jj].rs += (ii-2)*in.vs_drs;
      if (vs_in[jj].rs < 0.0) vs_in[jj].rs = 0.0;
    }
  }

  // State points with modified degeneracy parameter
  if (finite_temperature){
    for (int ii=0; ii<VSS_NUMEL; ii=ii+VSS_STENCIL){
      for(int jj=ii; jj<ii+VSS_STENCIL; jj++){
	vs_in[jj].Theta += (ii/VSS_STENCIL-2)*in.vs_dt;
	if (vs_in[jj].Theta < 0.0) vs_in[jj].Theta = 0.0;
      }
    }
  }

  // Chemical potential
  if (verbose) printf("Chemical potential calculation: ");
  if (finite_temperature)
    compute_chemical_potential_vs_stls(vs_in);
  if (verbose) printf("Done. \n");
  
  // Normalized ideal Lindhard density response
  if (verbose) printf("Normalized ideal Lindhard density calculation: ");
  if (finite_temperature)
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
  for (int ii=1; ii< in->vs_nrs; ii++) {
    rsa[ii] =  in->vs_drs*(ii);
  }
  
}

// ---------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE INITIAL GUESS
// ---------------------------------------------------------------------
void initial_guess_vs_stls(vs_struct xx, vs_struct SS, vs_struct SSHF,
			   vs_struct GG, vs_struct GG_new, vs_struct phi,
			   input *vs_in){

  input in = vs_in[VSS_IDXIN];
  
  if (strcmp(in.stls_guess_file,"NO_FILE")==0){

    for (int ii=0; ii<VSS_NUMEL; ii++){
      initial_guess_stls(xx.el[ii], SS.el[ii], SSHF.el[ii],
			 GG.el[ii], GG_new.el[ii], phi.el[ii],
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
double vs_stls_thermo_iterations(vs_struct xx, vs_struct rsu, vs_struct rsa,
				 input *vs_in, bool verbose) {

  input in = vs_in[VSS_IDXIN];
  double alpha = in.vs_alpha;
  double iter_err = 1.0;
  int iter_counter = 0;
  
  // Iterations
  while (iter_counter < in.nIter && iter_err > in.vs_err_min_iter ) {

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
      printf("alpha (CSR): %.5e\n", alpha);
      fflush(stdout);
    }
    
  }

  
  return alpha;
  
}

// Compute error for the iterations used to enforce the CSR
double vs_stls_thermo_err(double alpha, input *vs_in){

  return fabs((alpha - vs_in[VSS_IDXIN].vs_alpha)/alpha);
  
}

// Update self-consistency parameter for the CSR
void vs_stls_thermo_update(double alpha, input *vs_in){

  for (int ii=0; ii<VSS_NUMEL; ii++) {
    vs_in[ii].vs_alpha = alpha;
  }
   
}


// Iterations over the structural properties
void vs_stls_struct_iterations(vs_struct SS, vs_struct SSHF,
			       vs_struct GG, vs_struct GG_new,
			       vs_struct phi, vs_struct xx,
			       input *vs_in, bool verbose) {

  input in = vs_in[VSS_IDXIN];
    
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
  if (iter_err > in.err_min_iter) {
    fprintf(stderr, "The calculation for the structural properties for "
  	    "the state point (rs = %f, theta = %f) did not converge\n",
  	    vs_in[VSS_IDXIN].rs, vs_in[VSS_IDXIN].Theta);
    exit(EXIT_FAILURE);
  }
	    
  if (verbose) printf("Done.\n");
 
}

// Compute error for the iterations over the structural properties
double vs_stls_struct_err(vs_struct GG, vs_struct GG_new, input *vs_in){

  return stls_err(GG.rst, GG_new.rst, vs_in[VSS_IDXIN]);
  
}

// Update structural properties
void vs_stls_struct_update(vs_struct GG, vs_struct GG_new, input *vs_in){

  for (int ii=0; ii<VSS_NUMEL; ii++) {
    stls_update(GG.el[ii], GG_new.el[ii], vs_in[ii]);
  }
   
}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

void compute_chemical_potential_vs_stls(input *vs_in){

  #pragma omp parallel for
  for (int ii=0; ii<VSS_NUMEL; ii++){
    if (vs_in[ii].Theta > 0.0) {
      vs_in[ii].mu = compute_chemical_potential(vs_in[ii]);
    }
  }
  
}
					
// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_vs_stls(vs_struct SS, vs_struct SSHF, vs_struct GG, vs_struct phi,
		    vs_struct xx, input *vs_in){

  #pragma omp parallel for
  for (int ii=0; ii<VSS_NUMEL; ii++) {
    compute_ssf_stls(SS.el[ii], SSHF.el[ii], GG.el[ii],
		     phi.el[ii], xx.el[ii], vs_in[ii]);
  }
  
}

void compute_ssf_HF_vs_stls(vs_struct SS, vs_struct xx, input *vs_in){

  # pragma omp parallel for
  for (int ii=0; ii<VSS_NUMEL; ii++) {
    if (vs_in[ii].Theta == 0) {
      compute_ssf_HF_zero_temperature(SS.el[ii], xx.el[ii], vs_in[ii]);
    }
    else {
      compute_ssf_HF_finite_temperature(SS.el[ii], xx.el[ii], vs_in[ii]);
    }
	
  }
  
}
		    
// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_vs_stls(vs_struct GG, vs_struct SS, vs_struct xx, input *vs_in) {

  bool finite_temperature = false;
  vs_struct GG_stls;
  double *pp;
  input in = vs_in[VSS_IDXIN];

  // Check if finite temperature calculations must be performed
  if (in.Theta > 0.0) finite_temperature = true;
  
  // Allocate arrays
  for (int ii=0; ii<VSS_NUMEL; ii++){
    pp = malloc(sizeof(double) * in.nx);
    GG_stls.el[ii] = pp;
  }

  // STLS contribution to the VS-STLS static local field correction
  slfc_vs_stls_stls(GG, SS, GG_stls, xx, vs_in);

  // Wave-vector derivative contribution to the VS-STLS static local field correction
  slfc_vs_stls_dx(GG, GG_stls, xx, vs_in);
  
  // State point derivative contribution to the VS-STLS static local field correction
  slfc_vs_stls_drs(GG, GG_stls, vs_in);
  if (finite_temperature)
    slfc_vs_stls_dt(GG, GG_stls, vs_in);
  
  // Free memory
  for (int ii=0; ii<VSS_NUMEL; ii++){
    free(GG_stls.el[ii]);
  }
  
}

// STLS contribution to the VS-STLS static local field correction
void slfc_vs_stls_stls(vs_struct GG, vs_struct SS, vs_struct GG_stls,
		       vs_struct xx, input *vs_in){

  input in = vs_in[VSS_IDXIN];

  #pragma omp parallel for
  for (int ii=0; ii<VSS_NUMEL; ii++){
    compute_slfc(GG_stls.el[ii], SS.el[ii], xx.el[ii], vs_in[ii]);
    for (int jj=0; jj<in.nx; jj++){
      GG.el[ii][jj] = GG_stls.el[ii][jj];
    }
  }
  
}

// Wave-vector derivative contribution to the VS-STLS static local field correction
void slfc_vs_stls_dx(vs_struct GG, vs_struct GG_stls, vs_struct xx,
		     input *vs_in){

  input in = vs_in[VSS_IDXIN];
  double a_dx = in.vs_alpha/(3.0*2.0*in.dx);

  #pragma omp parallel for
  for (int ii=0; ii<VSS_NUMEL; ii++){
    for (int jj=1; jj<in.nx-1; jj++){
      GG.el[ii][jj] += -a_dx*xx.el[ii][jj]*(GG_stls.el[ii][jj+1]
					    - GG_stls.el[ii][jj-1]);
    }
  }
      
}


// State point derivateive contribution to the VS-STLS static local field correction
void slfc_vs_stls_drs(vs_struct GG, vs_struct GG_stls, input *vs_in){

  #pragma omp parallel for
  for (int ii=0; ii<VSS_NUMEL; ii++){
    for (int jj=0; jj<vs_in[ii].nx; jj++){

      // Forward difference for all state points with rs - 2*drs
      if (ii % VSS_STENCIL == 0)
	
	GG.el[ii][jj] += -(1.0/3.0)*vs_in[ii].vs_alpha
	  *slfc_vs_stls_drs_forward(GG_stls, ii, jj, vs_in);

      // Backward difference for all state points with rs + 2*drs
      else if ((ii-VSS_STENCIL+1) % VSS_STENCIL == 0)

	GG.el[ii][jj] += -(1.0/3.0)*vs_in[ii].vs_alpha
	  *slfc_vs_stls_drs_backward(GG_stls, ii, jj, vs_in);

      // Centered difference for all the other state points
      else

	GG.el[ii][jj] += -(1.0/3.0)*vs_in[ii].vs_alpha
	  *slfc_vs_stls_drs_centered(GG_stls, ii, jj, vs_in);
      
    }
  }
      
}

void slfc_vs_stls_dt(vs_struct GG, vs_struct GG_stls, input *vs_in){

  # pragma omp parallel for
  for (int ii=0; ii<VSS_NUMEL; ii++){
    for (int jj=0; jj<vs_in[ii].nx; jj++){

      // Forward difference for all state points with Theta - 2*dt
      if (ii/VSS_STENCIL == 0)
	
	GG.el[ii][jj] += -(2.0/3.0)*vs_in[ii].vs_alpha
	  *slfc_vs_stls_dt_forward(GG_stls, ii, jj, vs_in);

      // Backward difference for all state points with Theta + 2*dt
      else if (ii/VSS_STENCIL == VSS_STENCIL-1)

	GG.el[ii][jj] += -(2.0/3.0)*vs_in[ii].vs_alpha
	  *slfc_vs_stls_dt_backward(GG_stls, ii, jj, vs_in);

      // Centered difference for all the other state points
      else

	GG.el[ii][jj] += -(2.0/3.0)*vs_in[ii].vs_alpha
	  *slfc_vs_stls_dt_centered(GG_stls, ii, jj, vs_in);
      
    }
  }
   
      
}

double slfc_vs_stls_drs_centered(vs_struct GG, int ii, int jj, input *vs_in){

  return vs_in[ii].rs*(GG.el[ii+1][jj] - GG.el[ii-1][jj])/(2.0*vs_in[ii].vs_drs);
  
}

double slfc_vs_stls_drs_forward(vs_struct GG, int ii, int jj, input *vs_in){

  return vs_in[ii].rs*(-3.0*GG.el[ii][jj] + 4.0*GG.el[ii+1][jj]
		       - GG.el[ii+2][jj])/(2.0*vs_in[ii].vs_drs);
  
}

double slfc_vs_stls_drs_backward(vs_struct GG, int ii, int jj, input *vs_in){

  return vs_in[ii].rs*(3.0*GG.el[ii][jj] - 4.0*GG.el[ii-1][jj]
		       + GG.el[ii-2][jj])/(2.0*vs_in[ii].vs_drs);
  
}

double slfc_vs_stls_dt_centered(vs_struct GG, int ii, int jj, input *vs_in){

  return vs_in[ii].Theta*(GG.el[ii+VSS_STENCIL][jj]
			  - GG.el[ii-VSS_STENCIL][jj])/(2.0*vs_in[ii].vs_dt);
  
}

double slfc_vs_stls_dt_forward(vs_struct GG, int ii, int jj, input *vs_in){

  return vs_in[ii].Theta*(-3.0*GG.el[ii][jj] + 4.0*GG.el[ii+VSS_STENCIL][jj]
			  - GG.el[ii+2*VSS_STENCIL][jj])/(2.0*vs_in[ii].vs_dt);
  
}

double slfc_vs_stls_dt_backward(vs_struct GG, int ii, int jj, input *vs_in){

  return vs_in[ii].Theta*(3.0*GG.el[ii][jj] - 4.0*GG.el[ii-VSS_STENCIL][jj]
			  + GG.el[ii-2*VSS_STENCIL][jj])/(2.0*vs_in[ii].vs_dt);
  
}

// -----------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY RESPONSE
// -----------------------------------------------------------------------

void compute_idr_vs_stls(vs_struct phi, vs_struct xx, input *vs_in){

  #pragma omp parallel for 
  for (int ii=0; ii<VSS_NUMEL; ii++) {
    if (vs_in[ii].Theta > 0.0) {
      compute_idr(phi.el[ii], xx.el[ii], vs_in[ii], false);
    }
  }
  
}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE PARAMETER FOR THE COMPRESSIBILITY RULE
// -------------------------------------------------------------------

double compute_alpha(vs_struct xx, vs_struct rsu, vs_struct rsa, input *vs_in){

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
  input in = vs_in[VSS_IDXIN];
  int nrs = in.vs_nrs;
  
  // Check if finite temperature calculations must be performed
  if (in.Theta > 0.0) finite_temperature = true;
  
  // Get integrand for the free energy
  compute_rsu(xx, rsu, rsa, vs_in, false);

  // Internal energy
  ursm = rsu.rst[nrs-4]/rsa.rst[nrs-4];
  urs = rsu.rst[nrs-3]/rsa.rst[nrs-3];
  ursp = rsu.rst[nrs-2]/rsa.rst[nrs-2];
   if (finite_temperature) {
     utp = rsu.rstp1[nrs-3]/rsa.rstp1[nrs-3];
     utm = rsu.rstm1[nrs-3]/rsa.rstm1[nrs-3];
   } 

 
  // Free energy
  frs = compute_free_energy(rsu.rst, rsa.rst, vs_in[VSS_IDXIN]);
  frsp = compute_free_energy(rsu.rst, rsa.rst, vs_in[VSS_IDXIN+1]);
  frsm = compute_free_energy(rsu.rst, rsa.rst, vs_in[VSS_IDXIN-1]);
  if (finite_temperature) {
    ftp = compute_free_energy(rsu.rstp1, rsa.rstp1, vs_in[VSS_IDXIN]);
    ftm = compute_free_energy(rsu.rstm1, rsa.rstm1, vs_in[VSS_IDXIN]);
    frsptp = compute_free_energy(rsu.rstp1, rsa.rstp1, vs_in[VSS_IDXIN+1]);
    frsptm = compute_free_energy(rsu.rstm1, rsa.rstm1, vs_in[VSS_IDXIN+1]);
    frsmtp = compute_free_energy(rsu.rstp1, rsa.rstp1, vs_in[VSS_IDXIN-1]);
    frsmtm = compute_free_energy(rsu.rstm1, rsa.rstm1, vs_in[VSS_IDXIN-1]);
  }
  
  // Internal energy derivatives
  dudrs = (ursp - ursm)/(2.0*in.vs_drs);
  if (finite_temperature)
    dudt =  (utp - utm)/(2.0*in.vs_dt);
  
  // Free energy derivatives
  dfdrs = (frsp - frsm)/(2.0*in.vs_drs);
  d2fdrs2 = (frsp -2.0*frs + frsm)/(in.vs_drs*in.vs_drs);
  if (finite_temperature) {
    dfdt = (ftp - ftm)/(2.0*in.vs_dt);
    d2fdt2 = (ftp -2.0*frs + ftm)/(in.vs_dt*in.vs_dt);
    d2fdrsdt = (frsptp - frsmtp - frsptm + frsmtm)/(4.0*in.vs_drs*in.vs_dt);
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

void compute_rsu(vs_struct xx, vs_struct rsu, vs_struct rsa,
		 input *vs_in, bool verbose){

  int nrs = vs_in[VSS_IDXIN].vs_nrs;
  int cutoff_idx = 0;
  
  // Arrays for VS-STLS solution
  vs_struct xx_tmp;
  vs_struct phi;
  vs_struct GG;
  vs_struct GG_new;
  vs_struct SS;
  vs_struct SSHF;
  
  // Temporary arrays for allocation
  double *xx_el = NULL; 
  double *phi_el = NULL;
  double *GG_el = NULL;
  double *GG_new_el = NULL;
  double *SS_el = NULL;
  double *SSHF_el = NULL;
   
  // Allocate arrays
  for (int ii=0; ii<VSS_NUMEL; ii++) {
    alloc_stls_arrays(vs_in[VSS_IDXIN], &xx_el, &phi_el, &GG_el,
		      &GG_new_el, &SS_el, &SSHF_el);
    xx_tmp.el[ii] = xx_el;
    phi.el[ii] = phi_el;
    GG.el[ii] = GG_el;
    GG_new.el[ii] = GG_new_el;
    SS.el[ii] = SS_el;
    SSHF.el[ii] = SSHF_el;
  }

  // Read file for thermodynamic integration
  if (strcmp(vs_in[VSS_IDXIN].vs_thermo_file, "NO_FILE") != 0)
    read_thermo_vs_stls(rsa, rsu, &cutoff_idx, vs_in);

  // Update the first part of rsu in steps of five
  compute_rsu_blocks(SS, SSHF, GG, GG_new, phi, xx,
		     rsu, rsa, vs_in, &cutoff_idx,
		     cutoff_idx+2, nrs-2, VSS_STENCIL,
		     true, verbose);
 
  // Update the remainder of rsu in steps of one
  compute_rsu_blocks(SS, SSHF, GG, GG_new, phi, xx,
		     rsu, rsa, vs_in, &cutoff_idx,
		     cutoff_idx, nrs, 1,
		     false, verbose);

  // Final check to ensure that all the rsu array has been filled
  if (cutoff_idx < nrs-1) {
    fprintf(stderr, "Unexpected error, the free energy integrand was computed "
	    "only up to r_s = %f\n", rsa.rst[cutoff_idx]);
    exit(EXIT_FAILURE);
  }

  // Free memory
  for (int ii=0; ii<VSS_NUMEL; ii++){
    free_stls_arrays(xx_tmp.el[ii], phi.el[ii], GG.el[ii],
		     GG_new.el[ii], SS.el[ii], SSHF.el[ii]);
  }
  
}

void compute_rsu_blocks(vs_struct SS, vs_struct SSHF,
			vs_struct GG, vs_struct GG_new,
			vs_struct phi, vs_struct xx,
			vs_struct rsu, vs_struct rsa,
			input *vs_in, int *last,
			int start, int end, int step,
			bool compute_guess, bool verbose){

  int nrs = vs_in[VSS_IDXIN].vs_nrs;
  int last_tmp = -1;
  input vs_in_tmp[VSS_NUMEL];
  double u_int;

  for (int ii=start; ii<end; ii=ii+step) {
    
    // Define state point
    for (int jj=0; jj<VSS_NUMEL; jj++) {
      vs_in_tmp[jj] = vs_in[jj];
    }
    vs_in_tmp[VSS_IDXIN].rs = rsa.rst[ii];

    // Initialize arrays that depend only on the state point
    init_tmp_vs_stls_arrays(vs_in_tmp, xx, phi, SSHF, verbose);

    // Initial guess
    /* if (compute_guess && ii == start) { */
      initial_guess_vs_stls(xx, SS, SSHF, GG, GG_new, phi, vs_in_tmp);
    /* } */
    
    // Compute structural properties with the VS-STLS static local field correction
    vs_stls_struct_iterations(SS, SSHF, GG, GG_new, phi, xx, vs_in_tmp, verbose);
    
    // Compute the free energy integrand
    for (int jj=0; jj<VSS_NUMEL; jj++){
      
      // Avoid divergencies for rs = 0.0 (the internal energy scales as 1.0/rs)
      if (vs_in_tmp[jj].rs == 0.0) vs_in_tmp[jj].rs = 1.0;
      
      // Internal energy
      u_int = compute_internal_energy(SS.el[jj], xx.el[jj], vs_in_tmp[jj]);
      
      // Free energy integrand
      rsu.el[jj][ii] = vs_in_tmp[jj].rs*u_int;
      
    }

    // Update free energy integrand of neighboring points
    // We care to keep synchronized only the rst, rstm1 and rstp1 elements of rsu,
    // since these are the only elements used to compute alpha (see compute_alpha)
    if (start>=2 && end<=nrs-2){
      
      for (int jj=-VSS_STENCIL; jj<2*VSS_STENCIL; jj=jj+VSS_STENCIL){
	rsu.el[VSS_IDXIN+jj][ii-2] = rsu.el[VSS_IDXIN+jj-2][ii];
	rsu.el[VSS_IDXIN+jj][ii-1] = rsu.el[VSS_IDXIN+jj-1][ii];
	rsu.el[VSS_IDXIN+jj][ii+1] = rsu.el[VSS_IDXIN+jj+1][ii];
	rsu.el[VSS_IDXIN+jj][ii+2] = rsu.el[VSS_IDXIN+jj+2][ii];
      }

    }
    
    // Keep track of the index of the last element that was processed
    last_tmp = ii;

  }

  if (last_tmp >= 0) {
    if (start>=2 && end<=nrs-2)
      *last = last_tmp + 3;
    else
      *last = last_tmp + 1;
  }
  
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FREE ENERGY
// -------------------------------------------------------------------

struct fxc_params {

  gsl_spline *rsu_sp_ptr;
  gsl_interp_accel *rsu_acc_ptr;
  
};

double compute_free_energy(double *rsu, double *rsa, input in) {

  double err;
  size_t neval;
  double fre;

  // Declare accelerator and spline objects
  gsl_spline *rsu_sp_ptr;
  gsl_interp_accel *rsu_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  rsu_sp_ptr = gsl_spline_alloc(gsl_interp_linear, in.vs_nrs);
  rsu_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(rsu_sp_ptr, rsa, rsu, in.vs_nrs);

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &fxc;

  // Internal energy
  struct fxc_params fxcp = {rsu_sp_ptr, rsu_acc_ptr};
  ff_int.params = &fxcp;
  gsl_integration_cquad(&ff_int,
  			rsa[0], in.rs,
  			0.0, QUAD_REL_ERR,
  			wsp,
  			&fre, &err, &neval);
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(rsu_sp_ptr);
  gsl_interp_accel_free(rsu_acc_ptr);

  // Output
  return fre/(in.rs*in.rs);

}

double fxc(double rs, void* pp) {

  struct fxc_params* params = (struct fxc_params*)pp;
  gsl_spline* rsu_sp_ptr = (params->rsu_sp_ptr);
  gsl_interp_accel* rsu_acc_ptr = (params->rsu_acc_ptr);
  
  if (rs >= 0.0)
    return gsl_spline_eval(rsu_sp_ptr, rs, rsu_acc_ptr);
  else
    return 0.0;
  
}


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

// write text files for output
void write_text_vs_stls(double *SS, double *GG, double *phi,
			double *SSHF, double *xx, double *rsu,
			double *rsa, input in){

  // STLS arrays
  write_text_stls(SS, GG, phi, SSHF, xx, in);
  
  // Free energy output
  write_text_fxc(rsu, rsa, in);

  // Free parameter output
  write_text_alpha_CSR(in);
  
}


// Write free energy to text file
void write_text_fxc(double *rsu, double *rsa, input in){

  FILE* fid;
  char out_name[100];
  
  sprintf(out_name, "fxc_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the free energy\n");
    exit(EXIT_FAILURE);
  }
  
  fprintf(fid, "%.8e %.8e %.8e\n", in.rs, in.Theta, compute_free_energy(rsu, rsa, in));
  
  fclose(fid);

}


// Write free parameter for the compressibility sum rule to file
void write_text_alpha_CSR(input in){

    FILE* fid;
    char out_name[100];

    sprintf(out_name, "alpha_csr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        fprintf(stderr, "Error while creating the output file for the free parameter\n");
        exit(EXIT_FAILURE);
    }
    
    fprintf(fid, "%.8e %.8e %.8e\n", in.rs, in.Theta, in.vs_alpha);

    fclose(fid);

}

// Read guess from input file
void read_guess_vs_stls(vs_struct SS, vs_struct GG, input *vs_in){
  
  for (int ii=0;  ii<VSS_NUMEL; ii++) {
    read_guess_stls(SS.el[ii], GG.el[ii], vs_in[ii]);
  }
  
}

// write binary file to use for thermodynamic integration
void write_thermo_vs_stls(vs_struct rsa, vs_struct rsu, input *vs_in){

  // Input structure
  input in = vs_in[VSS_IDXIN];
  double Theta = in.Theta;
    
  // Name of output file
  char out_name[100];
  sprintf(out_name, "thermo_int_rs%.3f_theta%.3f_%s.bin", in.rs, in.Theta, in.theory);

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(out_name, "wb");
  if (fid == NULL) {
    fprintf(stderr,"Error while creating file for thermodynamic integration\n");
    exit(EXIT_FAILURE);
  }

  // Data for the grid used for thermodynamic integration
  fwrite(&in.vs_nrs, sizeof(int), 1, fid);
  fwrite(&in.vs_drs, sizeof(double), 1, fid);

  // Data for the degeneracy parameters
  fwrite(&Theta, sizeof(double), 1, fid);
  fwrite(&in.vs_dt, sizeof(double), 1, fid);
  
  // Grid used for thermodynamic integration
  fwrite(rsa.rstm1, sizeof(double), in.vs_nrs, fid);
  fwrite(rsa.rst, sizeof(double), in.vs_nrs, fid);
  fwrite(rsa.rstp1, sizeof(double), in.vs_nrs, fid);
  
  // Free energy integrand
  fwrite(rsu.rstm1, sizeof(double), in.vs_nrs, fid);
  fwrite(rsu.rst, sizeof(double), in.vs_nrs, fid);
  fwrite(rsu.rstp1, sizeof(double), in.vs_nrs, fid);
  
  // Close binary file
  fclose(fid);

}


// read binary file to use for thermodynamic integration
void read_thermo_vs_stls(vs_struct rsa, vs_struct rsu,
			 int *cutoff, input *vs_in){

  // Variables
  input in = vs_in[VSS_IDXIN];
  input in_file = vs_in[VSS_IDXIN];
  size_t it_read;
  int nrs_file;
  double drs_file;
  double dt_file;
  double Theta_file;
  double *rsu_el = NULL;
  double *rsa_el = NULL;
  vs_struct rsa_file;
  vs_struct rsu_file;
  
  // Open binary file
  FILE *fid = NULL;
  fid = fopen(in.vs_thermo_file, "rb");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file for thermodynamic integration\n");
    exit(EXIT_FAILURE);
  }

  // Initialize number of items read from input file
  it_read = 0;

  // Check that the data in the file for thermodynamic integration is consistent
  it_read += fread(&nrs_file, sizeof(int), 1, fid);
  it_read += fread(&drs_file, sizeof(double), 1, fid);
  it_read += fread(&Theta_file, sizeof(double), 1, fid);
  it_read += fread(&dt_file, sizeof(double), 1, fid);
  check_thermo_vs_stls(drs_file, dt_file, Theta_file,
		       in, it_read, 4, fid,
		       true, true, false);   

  // Allocate temporary arrays
  in_file.vs_nrs = nrs_file;
  for (int ii=0; ii<VSS_NUMEL; ii++) {
    alloc_vs_stls_arrays(in_file, &rsa_el, &rsu_el);
    rsa_file.el[ii] = rsa_el;
    rsu_file.el[ii] = rsu_el;
  }
  
  // Array of coupling parameters used for thermodynamic integration
  it_read += fread(rsa_file.rstm1, sizeof(double), nrs_file, fid);
  it_read += fread(rsa_file.rst, sizeof(double), nrs_file, fid);
  it_read += fread(rsa_file.rstp1, sizeof(double), nrs_file, fid);
  
  // Free energy integrand
  it_read += fread(rsu_file.rstm1, sizeof(double), nrs_file, fid);
  it_read += fread(rsu_file.rst, sizeof(double), nrs_file, fid);
  it_read += fread(rsu_file.rstp1, sizeof(double), nrs_file, fid);

  // Check that all items where read and the end-of-file was reached
  check_thermo_vs_stls(drs_file, dt_file, Theta_file,
		       in, it_read, 6*nrs_file + 4, fid,
		       false, true, true);   
 
  // Close binary file
  fclose(fid);

  // Output
  for (int ii=0; ii<nrs_file; ii++){

    if (ii>in.vs_nrs-1) 
      break;

    rsu.rstm1[ii] = rsu_file.rstm1[ii];
    rsu.rst[ii] = rsu_file.rst[ii];
    rsu.rstp1[ii] = rsu_file.rstp1[ii];
    
  }

  *cutoff = nrs_file;
  
  // Free memory
  for (int ii=0; ii<VSS_NUMEL; ii++){
    free_vs_stls_arrays(rsa_file.el[ii], rsu_file.el[ii]);
  }
 
}


// Check consistency of the thermodynamic integration data
void check_thermo_vs_stls(double drs, double dt, double Theta, input in,
			  size_t it_read, size_t it_expected, FILE *fid,
			  bool check_grid, bool check_items, bool check_eof){

  int buffer;
  double tol = 1e-10;
  
  // Check that the grid in the guess data is consistent with input
  if (check_grid) {
    
    if (fabs(drs-in.vs_drs) > tol || fabs(dt-in.vs_dt) > tol){
      fprintf(stderr,"Grid from thermodynamic integration file is incompatible with input\n");
      fprintf(stderr,"Resolution (drs)  : %.16f (input), %.16f (file)\n", in.vs_drs, drs);
      fprintf(stderr,"Resolution (dt)  : %.16f (input), %.16f (file)\n", in.vs_dt, dt);
      fclose(fid);
      exit(EXIT_FAILURE);
    }

    if (fabs(Theta-in.Theta) > tol){
      fprintf(stderr,"Degeneracy parameter from thermodynamic integration file is incompatible with input\n");
      fprintf(stderr,"Theta  : %.16f (input), %.16f (file)\n", in.Theta, Theta);
      fclose(fid);
      exit(EXIT_FAILURE);
    }

  }

  // Check that all the expected items where read
  if (check_items) {
    if (it_read != it_expected ) {
      fprintf(stderr,"Error while reading file for thermodynamic integration.\n");
      fprintf(stderr,"%ld Elements expected, %ld elements read\n", it_read, it_expected);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
  }
  
  // Check for end of file
  if (check_eof){
    it_read = fread(&buffer, sizeof(int), 1, fid); // Trigger end-of-file activation
    if (!feof(fid)) {
      fprintf(stderr,"Error while reading file for thermodynamic integration.\n");
      fprintf(stderr,"Expected end of file, but there is still data left to read.\n");
      fclose(fid);
      exit(EXIT_FAILURE);
    }
  }
  
  
}
