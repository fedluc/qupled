#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "stls.h"


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS EQUATIONS
// -------------------------------------------------------------------

void solve_stls(input in, bool verbose, 
		float **xx_out, float **SS_out, 
		float **SSHF_out, float **GG_out, 
		float **GG_new_out, float **phi_out) {

  // Arrays for STLS solution
  float *xx = NULL; 
  float *phi = NULL;
  float *GG = NULL;
  float *GG_new = NULL;
  float *SS = NULL;
  float *SSHF = NULL;

  // Start brand new calculation or load data from file
  bool init_flag = true;
  if (strcmp(in.phi_file,"NO_FILE") != 0) init_flag = false;
  if (init_flag) {
    alloc_stls_arrays(in, &xx, &phi, &GG, &GG_new, &SS, &SSHF);
  }
  else {
    read_bin(&in, &xx, &phi, &GG, &GG_new, &SS, &SSHF);
  }

  // Print on screen the parameter used to solve the STLS equation
  printf("------ Parameters used in the solution -------------\n");
  printf("Quantum degeneracy parameter: %f\n", in.Theta);
  printf("Quantum coupling parameter: %f\n", in.rs);
  printf("Chemical potential (low and high bound): %f %f\n",
	 in.mu_lo, in.mu_hi);
  printf("Wave-vector cutoff: %f\n", in.xmax);
  printf("Wave-vector resolutions: %f\n", in.dx);
  printf("Number of Matsubara frequencies: %d\n", in.nl);
  printf("Maximum number of iterations: %d\n", in.nIter);
  printf("Error for convergence: %.5e\n", in.err_min_iter);
  printf("----------------------------------------------------\n");
 

  if (init_flag){

    // Chemical potential
    if (verbose) printf("Chemical potential calculation: ");
    in.mu = compute_mu(in);
    if (verbose) printf("Done. Chemical potential: %.8f\n", in.mu);

    // Wave-vector grid
    if (verbose) printf("Wave-vector grid initialization: ");
    wave_vector_grid(xx, in);
    if (verbose) printf("Done.\n");
    
    // Normalized ideal Lindhard density
    if (verbose) printf("Normalized ideal Lindhard density calculation:\n");
    compute_phi(phi, xx, in, verbose);
    if (verbose) printf("Done.\n");
    
    // Static structure factor in the Hartree-Fock approximation
    if (verbose) printf("Static structure factor in the Hartree-Fock approximation: ");
    compute_ssfHF(SSHF, xx, in);
    if (verbose) printf("Done.\n");

  }

  // Initial guess for Static structure factor (SSF) and static-local field correction (SLFC)
  for (int ii=0; ii < in.nx; ii++) {
    GG[ii] = 0.0;
    GG_new[ii] = 1.0;
  }
  compute_ssf(SS, SSHF, GG, phi, xx, in);
  
  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  float iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    clock_t tic = clock();
    
    // Update SLFC
    compute_slfc(GG_new, SS, xx, in);
    
    // Update diagnostic
    iter_err = 0.0;
    iter_counter++;
    for (int ii=0; ii<in.nx; ii++) {
      iter_err += (GG_new[ii] - GG[ii]) * (GG_new[ii] - GG[ii]);
      GG[ii] = in.a_mix*GG_new[ii] + (1-in.a_mix)*GG[ii];
    }
    iter_err = sqrt(iter_err);
    
    // Update SSF
    compute_ssf(SS, SSHF, GG, phi, xx, in);
    
    // End timing
    clock_t toc = clock();
    
    // Print diagnostic
    if (verbose) {
      printf("--- iteration %d ---\n", iter_counter);
      printf("Elapsed time: %f seconds\n", ((float)toc - (float)tic) / CLOCKS_PER_SEC);
      printf("Residual error: %.5e\n", iter_err);
      fflush(stdout);
    }
  }
  if (verbose) printf("Done.\n");
  
  // Internal energy
  if (verbose) printf("Internal energy: %f\n",compute_uex(SS, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text(SS, GG, xx, in);
  if (init_flag) write_bin(phi, SSHF, in); 
  if (verbose) printf("Done.\n");

  // Output to variable
  bool free_xx = true, free_SS = true, free_SSHF = true, 
    free_GG = true, free_GG_new = true, free_phi = true;
  if (xx_out != NULL){ 
    *xx_out = xx;
    free_xx = false;
  }
  if (SS_out != NULL){
    *SS_out = SS;
    free_SS = false;
  }
  if (SSHF_out != NULL){
    *SSHF_out = SSHF;
    free_SSHF = false;
  }
  if (GG_out != NULL){
    *GG_out = GG;
    free_GG = false;
  }
  if (GG_new_out != NULL){
    *GG_new_out = GG_new;
    free_GG_new = false;
  }
  if (phi_out != NULL){
    *phi_out = phi;
    free_phi = false;
  }

  // Free memory
  free_stls_arrays(xx, free_xx, phi, free_phi, 
		   GG, free_GG, GG_new, free_GG_new,
		   SS, free_SS, SSHF, free_SSHF);

 
}

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays(input in, float **xx, float **phi, 
		       float **GG, float **GG_new, 
		       float **SS, float **SSHF){

  *xx = malloc( sizeof(float) * in.nx);
  *phi = malloc( sizeof(float) * in.nx * in.nl);
  *SSHF = malloc( sizeof(float) * in.nx);
  *GG = malloc( sizeof(float) * in.nx);
  *GG_new = malloc( sizeof(float) * in.nx);
  *SS = malloc( sizeof(float) * in.nx);
  
}

void free_stls_arrays(float *xx, bool free_xx, 
		      float *phi, bool free_phi, 
		      float *GG, bool free_GG, 
		      float *GG_new, bool free_GG_new, 
		      float *SS, bool free_SS,
		      float *SSHF, bool free_SSHF){

  if (free_xx) free(xx);
  if (free_phi) free(phi);
  if (free_SSHF) free(SSHF);
  if (free_SS) free(SS);
  if (free_GG) free(GG);
  if (free_GG_new) free(GG_new);
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

struct nc_params {

  float Theta;

};

float compute_mu(input in) {
  
  // Variables
  int max_iter = 100;
  double mu_lo = in.mu_lo;
  double mu_hi = in.mu_hi;
  double mu;
  int status, iter = 0;

  // Set-up function
  gsl_function ff_root;
  ff_root.function = &normalization_condition;
  struct nc_params ncp = {in.Theta};
  ff_root.params = &ncp;

  // Set-up root-solver
  const gsl_root_fsolver_type * rst = gsl_root_fsolver_bisection;
  gsl_root_fsolver * rs = gsl_root_fsolver_alloc(rst);
  gsl_root_fsolver_set(rs, &ff_root, mu_lo, mu_hi);

  // Solve normalization condition to find chemical potential
  do
  {
    
    // Solver iteration
    status = gsl_root_fsolver_iterate(rs);

    // Get solver status
    mu = gsl_root_fsolver_root(rs);
    mu_lo = gsl_root_fsolver_x_lower(rs);
    mu_hi = gsl_root_fsolver_x_upper(rs);
    status = gsl_root_test_interval (mu_lo, mu_hi,
				     0, 1e-10);
    
    // Update iteration counter
    iter++;

  }
  while (status == GSL_CONTINUE && iter < max_iter);

  // Free memory
  gsl_root_fsolver_free(rs);

  // Output
  return (float)mu;

}


double normalization_condition(double mu, void *pp) {

  struct nc_params *params = (struct nc_params*)pp;
  double Theta = (params->Theta);

  return gsl_sf_gamma(1.5)*gsl_sf_fermi_dirac_half(mu) 
    - 2.0/(3.0*pow(Theta, 3.0/2.0));

}

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void wave_vector_grid(float *xx, input in){
 
  xx[0] = in.dx/2.0;
  for (int ii=1; ii < in.nx; ii++) xx[ii] = xx[ii-1] + in.dx;

}

// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A TWO-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx2(int xx, int yy, int x_size) {
  return (yy * x_size) + xx;
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY
// -------------------------------------------------------------------

void compute_phi(float *phi, float *xx,  input in, bool verbose) {

  // Temporary array to store results
  float *phil = malloc( sizeof(float) * in.nx);
  
  // Loop over the Matsubara frequency
  for (int ll=0; ll<in.nl; ll++){
    if (verbose) printf("l = %d\n", ll);
    // Compute lindhard density
    compute_phil(phil, xx, ll, in);
    // Fill output array
    for (int ii=0; ii<in.nx; ii++){
      phi[idx2(ii,ll,in.nx)] = phil[ii];
    }    
  }
  
  // Free memory
  free(phil);

}

void compute_phil(float *phil, float *xx,  int ll, input in) {

  // Normalized ideal Lindhard density 
  for (int ii = 0; ii<in.nx; ii++) {

    phil[ii] = 0.0;
    if (ll == 0){
      
      for (int jj=0; jj<in.nx; jj++){
    	phil[ii] += phix0(xx[jj], xx[ii], in);
      }
      phil[ii] *= in.dx;
      
    }
    else {

      for (int jj=0; jj<in.nx; jj++){
    	phil[ii] += phixl(xx[jj], xx[ii], ll, in);
      }
      phil[ii] *= in.dx;
      

    }
    
  }
  
}

float phixl(float yy, float xx, int ll, input in) {

  float yy2 = yy*yy, xx2 = xx*xx, txy = 2*xx*yy, 
    tplT = 2*M_PI*ll*in.Theta, tplT2 = tplT*tplT;
  
  return 1.0/(2*xx)*yy/(exp(yy2/in.Theta - in.mu) + 1.0)
    *log(((xx2+txy)*(xx2+txy) + tplT2)/((xx2-txy)*(xx2-txy) + tplT2));


}

float phix0(float yy, float xx, input in) {

  float yy2 = yy*yy, xx2 = xx*xx, xy = xx*yy;

  if (xx < 2*yy){
    return 1.0/(in.Theta*xx)*((yy2 - xx2/4.0)*log((2*yy + xx)/(2*yy - xx)) + xy)
      *yy/(exp(yy2/in.Theta - in.mu) + exp(-yy2/in.Theta + in.mu) + 2.0);
  }
  else if (xx > 2*yy){
    return 1.0/(in.Theta*xx)*((yy2 - xx2/4.0)*log((2*yy + xx)/(xx - 2*yy)) + xy)
      *yy/(exp(yy2/in.Theta - in.mu) + exp(-yy2/in.Theta + in.mu) + 2.0);
  }
  else {
    return 1.0/(in.Theta)*yy2/(exp(yy2/in.Theta - in.mu) 
			       + exp(-yy2/in.Theta + in.mu) + 2.0);;
  }


}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

float ssfHF(float yy, float xx, input in) {

  float yy2 = yy*yy, ypx = yy + xx, ymx = yy - xx;
 
  return -3.0*in.Theta/(4.0*xx)*yy/(exp(yy2/in.Theta - in.mu) + 1.0)
    *log((1 + exp(in.mu - ymx*ymx/in.Theta))
	 /(1 + exp(in.mu - ypx*ypx/in.Theta)));

}

void compute_ssfHF(float *SS,  float *xx,  input in){

  // Static structure factor in the Hartree-Fock approximation
  for (int ii = 0; ii < in.nx; ii++) {

    SS[ii] = 0.0;
    for (int jj=0; jj<in.nx; jj++){
      SS[ii] += ssfHF(xx[jj], xx[ii], in);
    }
    SS[ii] *= in.dx;
    SS[ii] += 1.0;

  }
  
}


void compute_ssf(float *SS, float *SSHF, float *GG, 
		 float *phi, float *xx, input in){

  float lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  float pilambda = M_PI*lambda;
  float ff = 4*lambda*lambda*in.rs;
  float ff3_2T = 3.0*in.Theta*ff/2.0;
  float xx2, BB, BB_tmp, phixl;
  for (int ii=0; ii<in.nx; ii++){

    xx2 = xx[ii]*xx[ii];
    BB = 0.0;
    if (xx[ii] > 0.0){

      for (int ll=0; ll<in.nl; ll++){
	phixl = phi[idx2(ii,ll,in.nx)];
	BB_tmp = phixl*phixl/(pilambda*xx2 + ff*(1 - GG[ii])*phixl);
	if (ll>0) BB_tmp *= 2.0;
	BB += BB_tmp;
      }

      SS[ii] = SSHF[ii] - ff3_2T*(1- GG[ii])*BB;

    }
    else {
      SS[ii] = 0.0;
    }

  }

}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc(float *GG, float *SS, float *xx, input in) {

  // Static local field correction
  for (int ii = 0; ii < in.nx; ii++) {
    
    GG[ii] = 0.0;
    for (int jj=0; jj<in.nx; jj++){
      GG[ii] += slfc(xx[jj], xx[ii], SS[jj]);
    }
    GG[ii] *= in.dx;

  }
  
}

float slfc(float yy, float xx, float SS) {

    float yy2 = yy * yy, xx2 = xx * xx;

    if (xx > yy){
      return -3.0/4.0* yy2 * (SS - 1.0)
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(xx - yy)));
    }
    else if (xx < yy) {
      return -3.0/4.0* yy2 * (SS - 1.0)
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(yy - xx)));
    }
    else {
      return yy2 * (SS - 1.0);
    }

}



// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

float compute_uex(float *SS, input in) {

  float ie;
  float lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);  

  // Internal energy
  ie = 0.0;
  for (int jj=0; jj<in.nx; jj++){
    ie += SS[jj] - 1.0;
  }
  ie *= in.dx;
  
  // Output
  return ie/(M_PI*in.rs*lambda);

}

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------


// write text files with SSF and SLFC
void write_text(float *SS, float *GG, float *xx, input in){


    FILE* fid;
    
    // Output for SSF
    fid = fopen("ssf_STLS.dat", "w");
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
    fid = fopen("slfc_STLS.dat", "w");
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


// write binary file with density response
void write_bin(float *phi, float *SSHF, input in){


  // Open binary file
  FILE *fid = NULL;
  fid = fopen("dens_response.bin", "wb");
  if (fid == NULL) {
    fprintf(stderr,"Error while creating file with density response");
    exit(EXIT_FAILURE);
  }

  // Input data
  fwrite(&in, sizeof(input), 1, fid);

  // Density response
  fwrite(phi, sizeof(float), in.nx*in.nl, fid);

  // Static structure factor in the Hartree-Fock approximation
  fwrite(SSHF, sizeof(float), in.nx, fid);

  // Close binary file
  fclose(fid);

}
 

// read binary file with density response information
void read_bin(input *in, float **xx, float **phi, 
	      float **GG, float **GG_new, 
	      float **SS, float **SSHF){

  // Variables
  float *xx_local = NULL; 
  float *phi_local = NULL;
  float *GG_local = NULL;
  float *GG_new_local = NULL;
  float *SS_local = NULL;
  float *SSHF_local = NULL;
  input in_load;

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(in->phi_file, "rb");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file with density response");
    exit(EXIT_FAILURE);
  }

  // Input data from binary file
  fread(&in_load, sizeof(input), 1, fid);

  // Allocate arrays
  alloc_stls_arrays(in_load, &xx_local, &phi_local, 
		    &GG_local, &GG_new_local, &SS_local, 
		    &SSHF_local);

  // Chemical potential
  in_load.mu = compute_mu(in_load);

  // Wave-vector grid
  wave_vector_grid(xx_local, in_load);
  
  // Density response
  fread(phi_local, sizeof(float), in_load.nx*in_load.nl, fid);

  // Static structure factor in the Hartree-Fock approximation
  fread(SSHF_local, sizeof(float), in_load.nx, fid);

  // Close binary file
  fclose(fid);
  
  // Assign output
  in->Theta = in_load.Theta;
  in->dx = in_load.dx;
  in->xmax = in_load.xmax;
  in->nx = in_load.nx;
  in->nl = in_load.nl;  
  in->mu = in_load.mu;
  *xx = xx_local;
  *phi = phi_local;
  *GG = GG_local;
  *GG_new = GG_new_local;
  *SS = SS_local;
  *SSHF = SSHF_local;
	    
}
