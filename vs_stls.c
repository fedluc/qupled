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
  alloc_vs_stls_arrays(in, &uint, &rsArray);
  
  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_vs_stls_arrays(&in, xx, rsArray, verbose);
  init_state_point_vs_stls_arrays(&in, xx, phi, SSHF, verbose);
  
  // Initial guess
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

  // Iterative procedure to fix the compressibility sum rule
  if (verbose) printf("Iterative procedure to enforce the compressibility sum rule ...\n");
  in.a_csr = vs_stls_thermo_iterations(xx, uint, rsArray, in, true);
  if (verbose) printf("Done.\n");

  // Structural properties
  if (verbose) printf("Structural properties calculation...\n");
  vs_stls_struct_iterations(SS, SSHF, GG, GG_new, phi, xx, in, false);
  if (verbose) printf("Done.\n");
  
  // Thermodynamic properties
  if (verbose) printf("Free parameter: %.5f\n",in.a_csr);
  if (verbose) printf("Internal energy: %.10f\n",compute_internal_energy(SS, xx, in));
  if (verbose) printf("Free energy: %.10f\n",compute_free_energy(uint, rsArray, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text_stls(SS, GG, phi, SSHF, xx, in);
  write_text_vs_stls(uint, rsArray, in);
  write_guess_stls(SS, GG, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
  free_vs_stls_arrays(uint, rsArray);
 
}


// ---------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE VS-STLS SCHEME
// ---------------------------------------------------------------------

// Iterations over the parameter used to enforce the CSR rule
double vs_stls_thermo_iterations(double *xx, double *uint,
				 double *rsArray, input in,
				 bool verbose) {

  double alpha = in.a_csr;
  double iter_err = 1.0;
  int iter_counter = 0;

  // Iterations
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {

    // Start timing
    double tic = omp_get_wtime();

    // Get parameter to enforce the compressibility sum rule
    alpha = get_alpha(xx, uint, rsArray, in);

    // Update diagnostic
    iter_err = fabs((alpha - in.a_csr)/alpha);
    iter_counter++;
    in.a_csr = alpha;

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


// Iterations over the structural properties
void vs_stls_struct_iterations(double *SS, double *SSHF,
			       double *GG, double *GG_new,
			       double *phi, double *xx,
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
    compute_vs_slfc(GG_new, SS, xx, in);

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
      printf("--- iteration (structure) %d ---\n", iter_counter);
      printf("Elapsed time: %f seconds\n", toc - tic);
      printf("Residual error: %.5e\n", iter_err);
      fflush(stdout);
    }
    
  }
  if (verbose) printf("Done.\n");
 
 
}



// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_vs_stls_arrays(input in, double **uint, double **rsArray){
  *uint = malloc( sizeof(double) * in.nrs);
  *rsArray = malloc( sizeof(double) * in.nrs); 
}

void free_vs_stls_arrays(double *uint, double *rsArray){
  free(uint);
  free(rsArray);
}


// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

// Initialize arrays that do not depend on iterations and state points
void init_fixed_vs_stls_arrays(input *in, double *xx, double *rsArray,
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
  wave_vector_grid(xx, in);
  if (verbose) printf("Done.\n");

  // Array of coupling parameters
  if (verbose) printf("Coupling parameter grid initialization: ");
  rs_grid(rsArray, in);
  if (verbose) printf("Done.\n");
  
}


// Initialize arrays that do not depend on the iterations, but that
// are a function of the state point
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


// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rsArray, input *in){

  for (int ii=0; ii<in->nrs; ii++) {
    rsArray[ii] =  in->drs*(ii+1);
  }
  
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_vs_slfc(double *GG, double *SS, double *xx, input in) {

  bool finite_temperature = false;
  double *GGrs = NULL;
  double *GGrsp = NULL;
  double *GGrsm = NULL;
  double *GGtp = NULL;
  double *GGtm = NULL;
  double rs = in.rs;
  double Theta = in.Theta;
  double der_coeff = 2.0;
  
  // Check if finite temperature calculations must be performed
  if (Theta > 0.0) finite_temperature = true;
  
  // Allocate arrays
  GGrs = malloc( sizeof(double) * in.nx);
  GGrsp = malloc( sizeof(double) * in.nx);
  GGrsm = malloc( sizeof(double) * in.nx);
  if (finite_temperature) {
    GGtp = malloc( sizeof(double) * in.nx);
    GGtm = malloc( sizeof(double) * in.nx);
  }
  
  // STLS local field correction at the state point of interest
  compute_slfc(GGrs, SS, xx, in);
  
  // STLS local field correction for coupling parameter (rs) derivative
  in.rs = rs + in.drs;
  compute_slfc(GGrsp, SS, xx, in);
  in.rs = rs - in.drs;
  if (in.rs <= 0.0) {
    // Use a first order forward difference approximation for the derivative
    in.rs = rs;
    der_coeff = 1.0;
  }
  compute_slfc(GGrsm, SS, xx, in);
  in.rs = rs;
  
  // STLS local field correction for degeneracy parameter (Theta)  derivative
  // NOTE: For this derivative we use the same resolution used for the rs derivative
  if (finite_temperature) {
    if (in.drs >= Theta) {
      printf("Degeneracy parameter derivative cannot be computed, choose different vs-drs parameter\n");
      exit(EXIT_FAILURE);
    }
    in.Theta = Theta + in.drs;
    compute_slfc(GGtp, SS, xx, in);
    in.Theta = Theta - in.drs;
    compute_slfc(GGtm, SS, xx, in);
  }
  
  // VS-STLS static local field correction
  for (int ii=0; ii<in.nx; ii++){

    // STLS contribution
    GG[ii] = GGrs[ii];

    // State point derivative contribution
    GG[ii] += -in.a_csr/3.0*rs*(GGrsp[ii] - GGrsm[ii])/(der_coeff*in.drs);
    if (finite_temperature)
      GG[ii] += -in.a_csr*2.0/3.0*Theta*(GGtp[ii] - GGtm[ii])/(2.0*in.drs);
    
    // Wave-vector derivative contribution
    if (ii > 0 && ii < in.nx-1) {
      GG[ii] += -in.a_csr/3.0*xx[ii]*(GGrs[ii+1] - GGrs[ii-1])/(2.0*in.dx);
    }

    
  }
  
  // Free memory
  free(GGrs);
  free(GGrsp);
  free(GGrsm);
  free(GGtp);
  free(GGtm);
  
}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE PARAMETER FOR THE COMPRESSIBILITY RULE
// -------------------------------------------------------------------

double get_alpha(double *xx, double *uint,
		 double *rsArray, input in) {

  bool finite_temperature = false;
  double *uinttp = NULL;
  double *uinttm = NULL;  
  double frs, frsp, frsm;
  double ftp, ftm;
  double frsptp, frsmtp, frsptm, frsmtm;
  double dudrs, dudt;
  double dfdrs, dfdt;
  double d2fdrs2, d2fdt2, d2fdrsdt;
  double numer;
  double denom;
  double alpha;
  double rs = in.rs;
  double Theta = in.Theta;

  // Check if finite temperature calculations must be performed
  if (Theta > 0.0) finite_temperature = true;
  
  // Allocate temporary arrays for finite temperature calculations
  if (finite_temperature) {
    uinttp = malloc( sizeof(double) * in.nrs);
    uinttm = malloc( sizeof(double) * in.nrs);
  }
  
  // Get internal energy for thermodynamic integration
  fill_internal_energy_array(xx, uint, rsArray, in, false);
  if (finite_temperature) {
    in.Theta = Theta + in.drs;
    fill_internal_energy_array(xx, uinttp, rsArray, in, false);
    in.Theta = Theta - in.drs;
    fill_internal_energy_array(xx, uinttm, rsArray, in, false);
  }
  
  // Free energy
  frs = compute_free_energy(uint, rsArray, in);
  in.rs = rs + in.drs;
  frsp = compute_free_energy(uint, rsArray, in);
  in.rs = rs - in.drs;
  frsm = compute_free_energy(uint, rsArray, in);
  in.rs = rs;
  if (finite_temperature) {
    ftp = compute_free_energy(uinttp, rsArray, in);
    ftm = compute_free_energy(uinttm, rsArray, in);
    in.rs = rs + in.drs;
    frsptp = compute_free_energy(uinttp, rsArray, in);
    frsptm = compute_free_energy(uinttm, rsArray, in);
    in.rs = rs - in.drs;
    frsmtp = compute_free_energy(uinttp, rsArray, in);
    frsmtm = compute_free_energy(uinttm, rsArray, in);
    in.rs = rs;
  }
  
  // Internal energy derivatives
  dudrs = (uint[in.nrs-1] - uint[in.nrs-3])/(2.0*in.drs);
  if (finite_temperature) 
    dudt =  (uinttp[in.nrs-2] - uinttm[in.nrs-2])/(2.0*in.drs);
  
  // Free energy derivatives
  dfdrs = (frsp - frsm)/(2.0*in.drs);
  d2fdrs2 = (frsp -2.0*frs + frsm)/(in.drs*in.drs);
  if (finite_temperature) {
    dfdt = (ftp - ftm)/(2.0*in.drs);
    d2fdt2 = (ftp -2.0*frs + ftm)/(in.drs*in.drs);
    d2fdrsdt = (frsptp - frsmtp - frsptm + frsmtm)/(4.0*in.drs*in.drs);
  }
    
  // Parameter for the compressibility sum rule
  numer = (2.0*frs - (1.0/6.0)*rs*rs*d2fdrs2
	   + (4.0/3.0)*rs*dfdrs);
  denom = (uint[in.nrs-2] + (1.0/3.0)*rs*dudrs);
  if (finite_temperature) {
    numer += -(2.0/3.0)*Theta*Theta*d2fdt2
             -(2.0/3.0)*Theta*rs*d2fdrsdt
             +(1.0/3.0)*Theta*dfdt;
    denom += (2.0/3.0)*Theta*dudt;
  }
  alpha = numer/denom;

  // Free memory
  free(uinttp);
  free(uinttm);
  
  // Output
  return alpha;
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE INTERNAL ENERGY FOR THERMODYNAMIC
// INTEGRATION
// -------------------------------------------------------------------

void fill_internal_energy_array(double *xx, double *uint, double *rsArray,
				input in, bool verbose){


  #pragma omp parallel
  {
    #pragma omp for // Parallel calculations
    for (int ii=0; ii<in.nrs; ii++) {

      // Arrays for STLS solution
      double *xx_tmp = NULL; 
      double *phi = NULL;
      double *GG = NULL;
      double *GG_new = NULL;
      double *SS = NULL;
      double *SSHF = NULL;
      input in_tmp = in;

      // State point
      in_tmp.rs = rsArray[ii];
      
      // Allocate arrays
      alloc_stls_arrays(in_tmp, &xx_tmp, &phi, &GG, &GG_new, &SS, &SSHF);
      
      // Initialize arrays that depend only on the state point
      init_state_point_vs_stls_arrays(&in_tmp, xx, phi, SSHF, verbose);
      
      // Solve state point
      vs_stls_struct_iterations(SS, SSHF, GG, GG_new, phi, xx, in_tmp, verbose);
      
      // Internal energy
      uint[ii] = compute_internal_energy(SS, xx, in_tmp);
      
      // Free memory
      free_stls_arrays(xx_tmp, phi, GG, GG_new, SS, SSHF);
      
    }

 }
  
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FREE ENERGY
// -------------------------------------------------------------------

struct fex_params {

  gsl_spline *uint_sp_ptr;
  gsl_interp_accel *uint_acc_ptr;

};

double compute_free_energy(double *uint, double *rsArray, input in) {

  double err;
  size_t neval;
  double fre;

  // Declare accelerator and spline objects
  gsl_spline *uint_sp_ptr;
  gsl_interp_accel *uint_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  uint_sp_ptr = gsl_spline_alloc(gsl_interp_linear, in.nrs);
  uint_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(uint_sp_ptr, rsArray, uint, in.nrs);

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
  			rsArray[0]+1.0e-14, in.rs,
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

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

// write text files for output
void write_text_vs_stls(double *uint, double *rsArray, input in){

    FILE* fid;
    char out_name[100];
    
    // Output for the internal energy used for thermodynamic integration
    sprintf(out_name, "uint_thermoint_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static structure factor (HF)");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nrs; ii++)
        fprintf(fid, "%.8e %.8e\n", rsArray[ii], uint[ii]);

    fclose(fid);

    // Output for the free energy
    sprintf(out_name, "fxc_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the free energy\n");
        exit(EXIT_FAILURE);
    }
    fprintf(fid, "%.8e\n", compute_free_energy(uint, rsArray, in));
    fclose(fid);

    // Output for the free parameter
    sprintf(out_name, "alpha_csr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the free energy\n");
        exit(EXIT_FAILURE);
    }
    fprintf(fid, "%.8e\n", in.a_csr);
    fclose(fid);

}

