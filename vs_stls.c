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
// FUNCTION USED TO LOOP THE DATA STRUCTURES
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
    
  }

  printf("Invalid structure element for VS data structures\n"); 
  exit(EXIT_FAILURE);

}

// -------------------------------------------------------------------
// A NEW FUNCTION TO ALLOCATE THE STLS ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays_new(input in, stls_pointers *pp){

  pp->xx = malloc( sizeof(double) * in.nx);
  if (pp->xx == NULL) {
    fprintf(stderr, "Failed to allocate memory for the grid\n");
    exit(EXIT_FAILURE);
  }
  
  pp->phi = malloc( sizeof(double) * in.nx * in.nl);
    if (pp->phi == NULL) {
    fprintf(stderr, "Failed to allocate memory for the ideal density response\n");
    exit(EXIT_FAILURE);
  }

  pp->SS = malloc( sizeof(double) * in.nx);
    if (pp->SS == NULL) {
    fprintf(stderr, "Failed to allocate memory for the static structure factor\n");
    exit(EXIT_FAILURE);
  }
  
  pp->SSHF = malloc( sizeof(double) * in.nx);
  if (pp->SSHF == NULL) {
    fprintf(stderr, "Failed to allocate memory for the static structure factor"
  	    " in the Hartree-Fock approximation\n");
    exit(EXIT_FAILURE);
  }
  
  pp->GG = malloc( sizeof(double) * in.nx);
  if (pp->GG == NULL) {
    fprintf(stderr, "Failed to allocate memory for the static local field correction\n");
    exit(EXIT_FAILURE);
  }
  
  pp->GG_new = malloc( sizeof(double) * in.nx);
  if (pp->GG_new == NULL) {
    fprintf(stderr, "Failed to allocate memory for the static local field correction\n");
    exit(EXIT_FAILURE);
  }

  
}

void alloc_vs_stls_arrays(input in, stls_pointers pp, vs_sp *xx, vs_sp *phi,
			  vs_sp *SS, vs_sp *SSHF, vs_sp *GG, vs_sp *GG_new,
			  double **rsu, double **rs_arr, int el){

  bool alloc_failed = false;
  
  // Allocate arrays for thermodynamic integration
  if (el == 0) {
    
    *rsu = malloc( sizeof(double) * in.nrs);
    if (*rsu == NULL) {
      fprintf(stderr, "Failed to allocate memory for the free energy integrand\n");
      alloc_failed = true;
    }
   
    *rs_arr = malloc( sizeof(double) * in.nrs);
    if (*rsu == NULL) {
      fprintf(stderr, "Failed to allocate memory for the coupling parameter array\n");
      alloc_failed = true;
    }
    
  }

  // Allocate arrays for structural properties
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
    GG_new-> in= pp.GG_new;
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

  default:
    alloc_failed = true;

  }

  if (alloc_failed) {
    printf("The allocation of the VS-STLS arrays failed\n"); 
    exit(EXIT_FAILURE);
  }
 
}


void free_vs_stls_arrays(vs_sp xx, vs_sp phi, vs_sp SS, vs_sp SSHF,
			 vs_sp GG, vs_sp GG_new, double *rsu,
			 double *rs_arr){

  for (int ii=0; ii<VS_SP_EL; ii++) {
    free_stls_arrays(get_el(xx,ii), get_el(phi,ii), get_el(GG,ii),
		     get_el(GG_new, ii), get_el(SS,ii),
		     get_el(SSHF,ii));
  }
  free(rsu);
  free(rs_arr);
 
}

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
  double *rsu = NULL;
  double *rs_arr = NULL;
  input vs_in[VS_SP_EL];

  // Limit to zero temperature
  if (in.Theta > 0.0) {
    printf("The VS-STLS scheme is only implemented in the ground state\n");
    exit(EXIT_FAILURE);
  }
  // Allocate arrays
  stls_pointers stls_pp;
  for (int ii=0; ii<VS_SP_EL; ii++) {
    alloc_stls_arrays_new(in, &stls_pp);
    alloc_vs_stls_arrays(in, stls_pp, &xx, &phi,
			 &SS, &SSHF, &GG, &GG_new,
			 &rsu, &rs_arr, ii);
  }
  
  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_vs_stls_arrays(&in, vs_in, xx, rs_arr, verbose);
  init_state_point_vs_stls_arrays(vs_in, xx, phi, SSHF, verbose);
    
  // Initial guess
  if (strcmp(in.stls_guess_file,"NO_FILE")==0){
    // Static local field correction
    for (int ii=0; ii < in.nx; ii++) {
      // Static local field correction
      for (int jj=0; jj<VS_SP_EL; jj++) {
	get_el(GG,jj)[ii] = 0.0;
	get_el(GG_new,jj)[ii] = 0.0;
      }
    }
    // Static structure factor
    for (int ii=0;  ii<VS_SP_EL; ii++) {
      compute_ssf_stls(get_el(SS,ii), get_el(SSHF,ii),
		       get_el(GG,ii), get_el(phi,ii),
		       get_el(xx,ii), vs_in[ii]);
    }
  }
  else {
    // Read from file
    for (int ii=0;  ii<VS_SP_EL; ii++) {
      read_guess_stls(get_el(SS,ii), get_el(GG,ii),
		      vs_in[ii]);
    }
  }

  // Iterative procedure to fix the compressibility sum rule
  if (verbose) printf("Iterative procedure to enforce the compressibility sum rule ...\n");
  in.a_csr = vs_stls_thermo_iterations(xx, rsu, rs_arr, vs_in, true);
  if (verbose) printf("Done.\n");

  /* // Structural properties */
  /* if (verbose) printf("Structural properties calculation...\n"); */
  /* vs_stls_struct_iterations(SS, SSHF, GG, GG_new, phi, xx, in, false); */
  /* if (verbose) printf("Done.\n"); */
  
  /* // Thermodynamic properties */
  /* if (verbose) printf("Free parameter: %.5f\n",in.a_csr); */
  /* if (verbose) printf("Internal energy: %.10f\n",compute_internal_energy(SS, xx, in)); */
  /* if (verbose) printf("Free energy: %.10f\n",compute_free_energy(rsu, rs_arr, in)); */
  
  /* // Output to file */
  /* if (verbose) printf("Writing output files...\n"); */
  /* write_text_stls(SS, GG, phi, SSHF, xx, in); */
  /* write_text_vs_stls(rsu, rs_arr, in); */
  /* write_guess_stls(SS, GG, in); */
  /* if (verbose) printf("Done.\n"); */

  // Free memory
  free_vs_stls_arrays(xx, phi, SS, SSHF, GG, GG_new, rsu, rs_arr);
  
}


/* // ------------------------------------------------------------------- */
/* // FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS */
/* // ------------------------------------------------------------------- */

/* void alloc_vs_stls_arrays(input in, double **rsu, double **rs_arr){ */
/*   *rsu = malloc( sizeof(double) * in.nrs); */
/*   *rs_arr = malloc( sizeof(double) * in.nrs); */

/*   // Add error check */
/* } */

/* void free_vs_stls_arrays(double *rsu, double *rs_arr){ */
/*   free(rsu); */
/*   free(rs_arr); */
/* } */


// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

// Initialize arrays that do not depend on iterations and state points
void init_fixed_vs_stls_arrays(input *in, input *vs_in,
			       vs_sp xx, double *rs_arr,
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
  rs_grid(rs_arr, in);
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

    switch(ii) {
    case 1:
      vs_in[ii].rs += in.drs;
      break;
    case 2:
      vs_in[ii].rs -= in.drs;
      break;
    case 3:
      vs_in[ii].rs += 2.0*in.drs;
      break;
    case 4:
      vs_in[ii].rs -= 2.0*in.drs;
      break;
    }

    // Limit the state points to physical values
    if (vs_in[ii].rs < 0.0) vs_in[ii].rs = 0.0;
			       
  }
  
  // Chemical potential
  /* if (in->Theta > 0) { */
  /*   if (verbose) printf("Chemical potential calculation: "); */
  /*   in->mu = compute_chemical_potential(*in); */
  /*   if (verbose) printf("Done. Chemical potential: %.8f\n", in->mu); */
  /* } */
  
  // Normalized ideal Lindhard density response
  /* if (in->Theta > 0) { */
  /*   if (verbose) printf("Normalized ideal Lindhard density calculation:\n"); */
  /*   compute_idr(phi.in, xx.in, *in, verbose); */
  /*   compute_idr(phi.rsp1, xx.rsp1, *in, verbose); */
  /*   if (verbose) printf("Done.\n"); */
  /* } */
  
  // Static structure factor in the Hartree-Fock approximation
  if (verbose) printf("Static structure factor in the Hartree-Fock approximation: ");
  /* if (in->Theta == 0) { */
    for (int ii=0; ii<VS_SP_EL; ii++) {
      compute_ssf_HF_zero_temperature(get_el(SSHF,ii), get_el(xx,ii), vs_in[ii]);
    }
  /* } */
  /* else { */
  /*   for (int ii=0; ii<VS_SP_EL; ii++) { */
  /*     compute_ssf_HF_zero_temperature(get_el(SSHF,ii), get_el(xx,ii), vs_in[ii]); */
  /*   } */
  /* } */
  if (verbose) printf("Done.\n");

}


// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rs_arr, input *in){

  rs_arr[0] = 0.0;
  for (int ii=1; ii<in->nrs; ii++) {
    rs_arr[ii] =  in->drs*(ii);
  }
  
}



// ---------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE VS-STLS SCHEME
// ---------------------------------------------------------------------

// Iterations over the parameter used to enforce the CSR rule
double vs_stls_thermo_iterations(vs_sp xx, double *rsu,
				 double *rs_arr, input *vs_in,
				 bool verbose) {

  double alpha = vs_in[0].a_csr;
  int max_iter = vs_in[0].nIter;
  double min_err = vs_in[0].err_min_iter;

  // Iterations
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < max_iter && iter_err > min_err ) {

    // Start timing
    double tic = omp_get_wtime();

    // Get parameter to enforce the compressibility sum rule
    alpha = compute_alpha(xx, rsu, rs_arr, vs_in);

    // Update diagnostic
    iter_err = fabs((alpha - vs_in[0].a_csr)/alpha);
    iter_counter++;

    // Update alpha for all state points
    for (int ii=0; ii<VS_SP_EL; ii++) {
      vs_in[ii].a_csr = alpha;
    }

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
void vs_stls_struct_iterations(vs_sp SS, vs_sp SSHF,
			       vs_sp GG, vs_sp GG_new,
			       vs_sp phi, vs_sp xx,
			       input *vs_in, bool verbose) {

  double a_mix = vs_in[0].a_mix;
  int max_iter = vs_in[0].nIter;
  double min_err = vs_in[0].err_min_iter;
  double nx = vs_in[0].nx;
    
  // Iterations
  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < max_iter && iter_err > min_err ) {
    
    // Start timing
    double tic = omp_get_wtime();
    
    // Update SSF
    for (int ii=0; ii<VS_SP_EL; ii++) {
      compute_ssf_stls(get_el(SS,ii), get_el(SSHF,ii),
		       get_el(GG,ii), get_el(phi,ii),
		       get_el(xx,ii), vs_in[ii]);
    }
    
    // Update SLFC
    compute_vs_slfc(GG_new, SS, xx, vs_in);
    
    // Update diagnostic
    iter_err = 0.0;
    iter_counter++;
    for (int ii=0; ii<nx; ii++) {
      // Perhaps consider the minimum between all state points to
      // be solved simultaneously?
      iter_err += (get_el(GG_new,0)[ii] - get_el(GG,0)[ii])
	          * (get_el(GG_new,0)[ii] - get_el(GG,0)[ii]);
      for (int jj=0; jj<VS_SP_EL; jj++) {
	get_el(GG,ii)[jj] = a_mix*get_el(GG_new,ii)[jj]
	                    + (1 - a_mix)*get_el(GG,ii)[jj];
      }
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
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_vs_slfc(vs_sp GG, vs_sp SS, vs_sp xx, input *vs_in) {

  /* bool finite_temperature = false; */
  /* double *GGrs = NULL; */
  /* double *GGrs_arr = NULL; */
  /* double *GGrsm = NULL; */
  /* double *GGtp = NULL; */
  /* double *GGtm = NULL; */
  /* double drs = vs_in[0].rs; */
  /* double der_coeff = 2.0; */
  
  // Check if finite temperature calculations must be performed
  /* if (Theta > 0.0) finite_temperature = true; */
  
  // Allocate arrays
  /* GGrs = malloc( sizeof(double) * in.nx); */
  /* GGrs_arr = malloc( sizeof(double) * in.nx); */
  /* GGrsm = malloc( sizeof(double) * in.nx); */
  /* if (finite_temperature) { */
  /*   GGtp = malloc( sizeof(double) * in.nx); */
  /*   GGtm = malloc( sizeof(double) * in.nx); */
  /* } */




  
  // WE MUST STORE THE OLD RESULTS SOMEWHERE



  
  /* // STLS contribution at all the state points that have to be updated simultaneously */
  /* for (int ii=0; ii<VS_SP_EL; ii++) { */
  /*   compute_slfc(get_el(GG,ii), get_el(SS,ii), get_el(xx, ii), vs_in[ii]); */
  /* } */

  /* // State point derivative contributions */
  /* for (int ii=0; ii<in.nx; ii++) { */
  /*   for (int jj=0; jj<VS_SP_EL; jj++) { */

  /*     switch(ii) { */
	
  /*     case 0: */
  /* 	get_el(GG,jj)[ii] += -in.a_csr/3.0*rs*(get_el( - GG.in[ii])/(2.0*in.drs); */
  /*     case 1: */
  /* 	return sp.rsp1; */
  /*     case 2: */
  /* 	return sp.rsm1; */
  /*     case 3: */
  /* 	return sp.rsp2; */
  /*     case 4: */
  /* 	return sp.rsm2; */
	
  /*     } */
  /*   } */
  /* } */

  /* // Add derivative contributions */
  /* for (int ii=0; ii<in.nx; ii++){ */

  /*   // State point derivative contribution */
  /*   GG.in[ii] += -in.a_csr/3.0*rs*(GG.rsp1[ii] - GG.in[ii])/(in.drs); */
  /*   GG.rsp1[ii] += -in.a_csr/3.0*(rs + in.drs)*(GG.rsp1[ii] - GG.in[ii])/(in.drs); */
  /*   /\* if (finite_temperature) *\/ */
  /*   /\*   GG[ii] += -in.a_csr*2.0/3.0*Theta*(GGtp[ii] - GGtm[ii])/(2.0*in.drs); *\/ */
    
  /*   // Wave-vector derivative contribution */
  /*   if (ii > 0 && ii < in.nx-1) { */
  /*     GG.in[ii] += -in.a_csr/3.0*xx.in[ii]*(GG.in[ii+1] - GG.in[ii-1])/(2.0*in.dx); */
  /*     GG.rps1[ii] += -in.a_csr/3.0*xx.rps1[ii]*(GG.rsp1[ii+1] - GG.rsp1[ii-1])/(2.0*in.dx); */
  /*   } */

    
  /* } */
  
}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE PARAMETER FOR THE COMPRESSIBILITY RULE
// -------------------------------------------------------------------

double compute_alpha(vs_sp xx, double *rsu,
		     double *rs_arr, input *vs_in) {

  /* bool finite_temperature = false; */
  /* double *rsutp = NULL; */
  /* double *rsutm = NULL; */
  double urs, ursp, ursm;
  /* double utp, utm; */
  double frs, frsp, frsm;
  /* double ftp, ftm; */
  /* double frsptp, frsmtp, frsptm, frsmtm; */
  /* double dudrs, dudt; */
  /* double dfdrs, dfdt; */
  /* double d2fdrs2, d2fdt2, d2fdrsdt; */
  double dudrs;
  double dfdrs;
  double d2fdrs2;
  double numer;
  double denom;
  double alpha;
  double rs = vs_in[0].rs;
  /* double Theta = vs_in[0].Theta; */
  double drs = vs_in[0].drs;
  int nrs = vs_in[0].nrs;
  
  // Check if finite temperature calculations must be performed
  //if (Theta > 0.0) finite_temperature = true;
  
  // Allocate temporary arrays for finite temperature calculations
  /* if (finite_temperature) { */
  /*   rsutp = malloc( sizeof(double) * in.nrs); */
  /*   rsutm = malloc( sizeof(double) * in.nrs); */
  /* } */
  
  // Get integrand for the free energy
  compute_rsu(xx, rsu, rs_arr, vs_in, false);
  /* if (finite_temperature) { */
  /*   in.Theta = Theta + in.drs; */
  /*   compute_rsu(xx, rsutp, rs_arr, in, false); */
  /*   in.Theta = Theta - in.drs; */
  /*   compute_rsu(xx, rsutm, rs_arr, in, false); */
  /* } */

  // Internal energy
  urs = rsu[nrs-2]/rs_arr[nrs-2];
  ursp = rsu[nrs-1]/rs_arr[nrs-1];
  ursm = rsu[nrs-3]/rs_arr[nrs-3];
  /* if (finite_temperature) { */
  /*   utp = rsutp[in.nrs-2]/rs; */
  /*   utm = rsutm[in.nrs-2]/rs; */
  /* } */
  
  // Free energy
  frs = compute_free_energy(rsu, rs_arr, vs_in[0]);
  frsp = compute_free_energy(rsu, rs_arr, vs_in[1]);
  frsm = compute_free_energy(rsu, rs_arr, vs_in[2]);
  /* in.rs = rs; */
  /* if (finite_temperature) { */
  /*   ftp = compute_free_energy(rsutp, rs_arr, in); */
  /*   ftm = compute_free_energy(rsutm, rs_arr, in); */
  /*   in.rs = rs + in.drs; */
  /*   frsptp = compute_free_energy(rsutp, rs_arr, in); */
  /*   frsptm = compute_free_energy(rsutm, rs_arr, in); */
  /*   in.rs = rs - in.drs; */
  /*   frsmtp = compute_free_energy(rsutp, rs_arr, in); */
  /*   frsmtm = compute_free_energy(rsutm, rs_arr, in); */
  /*   in.rs = rs; */
  /* } */
  
  // Internal energy derivatives
  dudrs = (ursp - ursm)/(2.0*drs);
  /* if (finite_temperature) */
  /*   dudt =  (utp - utm)/(2.0*in.drs); */
  // Internal energy derivatives
  /* dudrs = (rsu[in.nrs-1] - rsu[in.nrs-3])/(2.0*in.drs); */
  /* if (finite_temperature)  */
  /*   dudt =  (rsutp[in.nrs-2] - rsutm[in.nrs-2])/(2.0*in.drs); */
  
  // Free energy derivatives
  dfdrs = (frsp - frsm)/(2.0*drs);
  d2fdrs2 = (frsp -2.0*frs + frsm)/(drs*drs);
  /* if (finite_temperature) { */
  /*   dfdt = (ftp - ftm)/(2.0*in.drs); */
  /*   d2fdt2 = (ftp -2.0*frs + ftm)/(in.drs*in.drs); */
  /*   d2fdrsdt = (frsptp - frsmtp - frsptm + frsmtm)/(4.0*in.drs*in.drs); */
  /* } */
    
  // Parameter for the compressibility sum rule
  numer = (2.0*frs - (1.0/6.0)*rs*rs*d2fdrs2
	   + (4.0/3.0)*rs*dfdrs);
  denom = (urs + (1.0/3.0)*rs*dudrs);
  /* if (finite_temperature) { */
  /*   numer += -(2.0/3.0)*Theta*Theta*d2fdt2 */
  /*            -(2.0/3.0)*Theta*rs*d2fdrsdt */
  /*            +(1.0/3.0)*Theta*dfdt; */
  /*   denom += (2.0/3.0)*Theta*dudt; */
  /* } */
  alpha = numer/denom;

  // Free memory
  /* free(rsutp); */
  /* free(rsutm); */
  
  // Output
  return alpha;
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE INTEGRAND FOR THE FREE ENERGY
// -------------------------------------------------------------------

void compute_rsu(vs_sp xx, double *rsu, double *rs_arr,
		 input *vs_in, bool verbose){


  int nrs = vs_in[0].nrs;
  
  #pragma omp parallel
  {
    #pragma omp for // Parallel calculations
    for (int ii=0; ii<nrs; ii++) {

      // Arrays for VS-STLS solution
      vs_sp tmp1;
      vs_sp phi;
      vs_sp GG;
      vs_sp GG_new;
      vs_sp SS;
      vs_sp SSHF;
      double *tmp2 = NULL;
      double *tmp3 = NULL;
      input vs_in_tmp[VS_SP_EL];
      stls_pointers stls_pp;

      // Define state point
      for (int jj=0; jj<VS_SP_EL; jj++) {
	vs_in_tmp[jj] = vs_in[jj];
      }
      vs_in_tmp[0].rs = rs_arr[ii];
      
      // Allocate arrays
      for (int ii=0; ii<VS_SP_EL; ii++) {
	alloc_stls_arrays_new(vs_in[0], &stls_pp);
	alloc_vs_stls_arrays(vs_in[0], stls_pp, &tmp1, &phi,
			     &SS, &SSHF, &GG, &GG_new,
			     &tmp2, &tmp3, ii);
      }

      
      // Initialize arrays that depend only on the state point
      init_state_point_vs_stls_arrays(vs_in_tmp, xx, phi, SSHF, verbose);
 
      if (rs_arr[ii] == 0.0) {

	// Get the integrand from the HF static structure factor
	vs_in_tmp[0].rs = 1.0;
	rsu[ii] = compute_internal_energy(get_el(SSHF,0), get_el(xx,0), vs_in_tmp[0]);
	
      }
      else {
	
	// Compute structural properties with the VS-STLS static local field correction
	vs_stls_struct_iterations(SS, SSHF, GG, GG_new, phi, xx, vs_in_tmp, verbose);
      
	// Compute the integrand for the free energy from the static structure factor
	rsu[ii] = rs_arr[ii]*compute_internal_energy(get_el(SSHF,0), get_el(xx,0),
						     vs_in_tmp[0]);
	
      }
	
      // Free memory
      free_vs_stls_arrays(tmp1, phi, SS, SSHF, GG, GG_new, tmp2, tmp3);
      
    }

 }
  
}



/* // ------------------------------------------------------------------- */
/* // FUNCTIONS FOR OUTPUT AND INPUT */
/* // ------------------------------------------------------------------- */

/* // write text files for output */
/* void write_text_vs_stls(double *rsu, double *rs_arr, input in){ */

/*     FILE* fid; */
/*     char out_name[100]; */
    
/*     // Output for the internal energy used for thermodynamic integration */
/*     sprintf(out_name, "rsu_thermoint_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory); */
/*     fid = fopen(out_name, "w"); */
/*     if (fid == NULL) { */
/*         fprintf(stderr, "Error while creating the output file for the static structure factor (HF)"); */
/*         exit(EXIT_FAILURE); */
/*     } */
/*     for (int ii = 0; ii < in.nrs; ii++) */
/*         fprintf(fid, "%.8e %.8e\n", rs_arr[ii], rsu[ii]); */

/*     fclose(fid); */

/*     // Output for the free energy */
/*     sprintf(out_name, "fxc_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory); */
/*     fid = fopen(out_name, "w"); */
/*     if (fid == NULL) { */
/*         fprintf(stderr, "Error while creating the output file for the free energy\n"); */
/*         exit(EXIT_FAILURE); */
/*     } */
/*     fprintf(fid, "%.8e\n", compute_free_energy(rsu, rs_arr, in)); */
/*     fclose(fid); */

/*     // Output for the free parameter */
/*     sprintf(out_name, "alpha_csr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory); */
/*     fid = fopen(out_name, "w"); */
/*     if (fid == NULL) { */
/*         fprintf(stderr, "Error while creating the output file for the free parameter\n"); */
/*         exit(EXIT_FAILURE); */
/*     } */
/*     fprintf(fid, "%.8e\n", in.a_csr); */
/*     fclose(fid); */

/* } */

