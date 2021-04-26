#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "solvers.h"
#include "chemical_potential.h"
#include "stls.h"


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS EQUATIONS
// -------------------------------------------------------------------

void solve_stls(input in, bool verbose) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *GG_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;

  // Allocate arrays
  alloc_stls_arrays(in, &xx, &phi, &GG, &GG_new, &SS, &SSHF);

  // Initialize arrays that are not modified with the iterative procedure
  init_fixed_stls_arrays(&in, xx, phi, SSHF, verbose);
  
  // Initial guess for Static structure factor (SSF) and static-local field correction (SLFC)
  if (strcmp(in.guess_file,"NO_FILE")==0){
    for (int ii=0; ii < in.nx; ii++) {
      GG[ii] = 0.0;
      GG_new[ii] = 1.0;
    }
    compute_ssf(SS, SSHF, GG, phi, xx, in);
  }
  else {
    read_guess(SS, GG, in);
  }


  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    double tic = omp_get_wtime();
    
    // Update SSF
    compute_ssf(SS, SSHF, GG, phi, xx, in);

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
  if (verbose) printf("Done.\n");

  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);

 
}

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays(input in, double **xx, double **phi, 
		       double **GG, double **GG_new, 
		       double **SS, double **SSHF){

  *xx = malloc( sizeof(double) * in.nx);
  *phi = malloc( sizeof(double) * in.nx * in.nl);
  *SSHF = malloc( sizeof(double) * in.nx);
  *GG = malloc( sizeof(double) * in.nx);
  *GG_new = malloc( sizeof(double) * in.nx);
  *SS = malloc( sizeof(double) * in.nx);
  
}

void free_stls_arrays(double *xx, double *phi, double *GG, 
		      double *GG_new, double *SS,
		      double *SSHF){

  free(xx);
  free(phi);
  free(SSHF);
  free(SS);
  free(GG);
  free(GG_new);
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_stls_arrays(input *in, double *xx, 
			    double *phi, double *SSHF, bool verbose){

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
 
  // Chemical potential
  if (verbose) printf("Chemical potential calculation: ");
  in->mu = compute_mu(*in);
  if (verbose) printf("Done. Chemical potential: %.8f\n", in->mu);
  
  // Wave-vector grid
  if (verbose) printf("Wave-vector grid initialization: ");
  wave_vector_grid(xx, *in);
  if (verbose) printf("Done.\n");
  
  // Normalized ideal Lindhard density
  if (verbose) printf("Normalized ideal Lindhard density calculation:\n");
  compute_phi(phi, xx, *in, verbose);
  if (verbose) printf("Done.\n");
  
  // Static structure factor in the Hartree-Fock approximation
  if (verbose) printf("Static structure factor in the Hartree-Fock approximation: ");
  compute_ssfHF(SSHF, xx, *in);
  if (verbose) printf("Done.\n");

}

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void wave_vector_grid(double *xx, input in){
 
  xx[0] = 0.0;
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

struct phixl_params {

  double xx;
  double mu;
  double Theta;
  double ll;

};

void compute_phi(double *phi, double *xx,  input in, bool verbose) {

  // Temporary array to store results
  double *phil = malloc( sizeof(double) * in.nx);
  
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

void compute_phil(double *phil, double *xx,  int ll, input in) {

  double err;
  size_t nevals;

  // Integration workspace 
  gsl_integration_cquad_workspace *wsp 
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  if (ll == 0) ff_int.function = &phix0;
  else ff_int.function = &phixl;

  // Normalized ideal Lindhard density 
  for (int ii = 0; ii < in.nx; ii++) {
    
    struct phixl_params phixlp = {xx[ii], in.mu, in.Theta, ll};
    ff_int.params = &phixlp;
    gsl_integration_cquad(&ff_int, 
			  xx[0], xx[in.nx-1], 
			  0.0, 1e-5, 
			  wsp, 
			  &phil[ii], &err, &nevals);

  }

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

double phixl(double yy, void *pp) {

  struct phixl_params *params = (struct phixl_params*)pp;
  double xx = (params->xx);
  double mu = (params->mu);
  double Theta = (params->Theta);
  double ll = (params->ll);
  double yy2 = yy*yy, xx2 = xx*xx, txy = 2*xx*yy, 
    tplT = 2*M_PI*ll*Theta, tplT2 = tplT*tplT;
  
  if (xx > 0.0) {
    return 1.0/(2*xx)*yy/(exp(yy2/Theta - mu) + 1.0)
      *log(((xx2+txy)*(xx2+txy) + tplT2)/((xx2-txy)*(xx2-txy) + tplT2));
  }
  else {
    return 0;
  }

}

double phix0(double yy, void *pp) {

  struct phixl_params *params = (struct phixl_params*)pp;
  double xx = (params->xx);
  double mu = (params->mu);
  double Theta = (params->Theta);
  double yy2 = yy*yy, xx2 = xx*xx, xy = xx*yy;

  if (xx > 0.0){

    if (xx < 2*yy){
      return 1.0/(Theta*xx)*((yy2 - xx2/4.0)*log((2*yy + xx)/(2*yy - xx)) + xy)
        *yy/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0);
    }
    else if (xx > 2*yy){
      return 1.0/(Theta*xx)*((yy2 - xx2/4.0)*log((2*yy + xx)/(xx - 2*yy)) + xy)
        *yy/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0);
    }
    else {
      return 1.0/(Theta)*yy2/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0);;
    }
  }

  else{
    return (2.0/Theta)*yy2/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0);
  }

}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf(double *SS, double *SSHF, double *GG, 
		 double *phi, double *xx, input in){

  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double ff = 4*lambda*in.rs/M_PI;
  double xx2, BB, BB_tmp, BB_den, phixl;
  //double Axl, tplT;

  for (int ii=0; ii<in.nx; ii++){

    if (xx[ii] > 0.0){

      xx2 = xx[ii]*xx[ii];
      BB = 0.0;
      
      for (int ll=0; ll<in.nl; ll++){
	//tplT = 2*M_PI*ll*in.Theta;
	phixl = phi[idx2(ii,ll,in.nx)];
	//Axl = (4.0/3.0)*xx2/(tplT*tplT + xx2*xx2);
	BB_den = 1.0 + ff/xx2*(1 - GG[ii])*phixl;
	//BB_tmp = phixl*phixl/BB_den - Axl*Axl;
	BB_tmp = phixl*phixl/BB_den;
	if (ll>0) BB_tmp *= 2.0;
	BB += BB_tmp;
	
      }
      
      /* SS[ii] = SSHF[ii] */
      /*   - 3.0/2.0*ff/xx2*in.Theta*(1- GG[ii])*BB */
      /*   - 1.0/3.0*ff/xx2/in.Theta*(1 - GG[ii])* */
      /*   (1.0/sinh(xx2/(2*in.Theta))* */
      /*    1.0/sinh(xx2/(2*in.Theta)) + */
      /*    2.0*in.Theta/xx2* */
      /*    1.0/tanh(xx2/(2*in.Theta))); */
      SS[ii] = SSHF[ii]
	- 3.0/2.0*ff/xx2*in.Theta*(1- GG[ii])*BB;

    }
    else
      SS[ii] = 0.0;

  }

}

struct ssfHF_params {

  double xx;
  double mu;
  double Theta;

};

void compute_ssfHF(double *SS,  double *xx,  input in){

  double err;
  size_t nevals;

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &ssfHF;

  // Static structure factor in the Hartree-Fock approximation
  for (int ii = 0; ii < in.nx; ii++) {


    struct ssfHF_params ssfHFp = {xx[ii], in.mu, in.Theta};
    ff_int.params = &ssfHFp;
    gsl_integration_cquad(&ff_int,
			  xx[0], xx[in.nx-1],
			  0.0, 1e-5,
			  wsp,
			  &SS[ii], &err, &nevals);

    SS[ii] += 1.0;

  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
}

double ssfHF(double yy, void* pp) {

  struct phixl_params *params = (struct phixl_params*)pp;
  double xx = (params->xx);
  double mu = (params->mu);
  double Theta = (params->Theta);
  double yy2 = yy*yy, ypx = yy + xx, ymx = yy - xx;
 
  if (xx > 0.0){
    return -3.0*Theta/(4.0*xx)*yy/(exp(yy2/Theta - mu) + 1.0)
      *log((1 + exp(mu - ymx*ymx/Theta))/(1 + exp(mu - ypx*ypx/Theta)));
  }
  else {
    return -3.0/2.0*yy2/(1.0 + cosh(yy2/Theta - mu));
  }

}

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

struct slfc_params {

  double xx;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};


void compute_slfc(double *GG, double *SS, double *xx, input in) {

  double err;
  size_t nevals;
  
  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);    

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &slfc;

  // Static local field correction
  for (int ii = 0; ii < in.nx; ii++) {
    struct slfc_params slfcp = {xx[ii], ssf_sp_ptr, ssf_acc_ptr};
    ff_int.params = &slfcp;
    gsl_integration_cquad(&ff_int,
			  xx[0], xx[in.nx-1],
			  0.0, 1e-5,
			  wsp,
			  &GG[ii], &err, &nevals);
    GG[ii] *= -3.0/4.0;
  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(ssf_sp_ptr);
  gsl_interp_accel_free(ssf_acc_ptr);
  
}

double slfc(double yy, void* pp) {

  struct slfc_params* params = (struct slfc_params*)pp;
  double xx = (params->xx);
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);
  double yy2 = yy * yy, xx2 = xx * xx;

  if (xx > 0.0 && yy > 0.0){

    if (xx > yy){
      return -(3.0/4.0)* yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0)
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(xx - yy)));
    }
    else if (xx < yy) {
      return -(3.0/4.0) * yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0)
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(yy - xx)));
    }
    else {
      return -(3.0/4.0) * yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0);
    }

  }
  else
    return 0;



}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

struct uex_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};

double compute_uex(double *SS, double *xx,  input in) {

  double err;
  size_t neval;
  double ie;
  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);  
  
  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_linear, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &uex;

  // Internal energy
  struct uex_params uexp = {ssf_sp_ptr, ssf_acc_ptr};
  ff_int.params = &uexp;  
  gsl_integration_cquad(&ff_int,
			xx[0], xx[in.nx-1],
			0.0, 1e-5,
			wsp,
			&ie, &err, &neval);
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(ssf_sp_ptr);
  gsl_interp_accel_free(ssf_acc_ptr);

  // Output
  return ie/(M_PI*in.rs*lambda);

}

double uex(double yy, void* pp) {

  struct uex_params* params = (struct uex_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);

  return gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1;
}


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------


// write text files for output
void write_text(double *SS, double *GG, double *phi, 
		double *SSHF, double *xx, input in){


    FILE* fid;
    
    // Output for SSF
    char out_name[100];
    sprintf(out_name, "ssf_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static structure factor");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nx; ii++)
        fprintf(fid, "%.8e %.8e\n", xx[ii], SS[ii]);

    fclose(fid);

    // Output for SLFC
    sprintf(out_name, "slfc_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static local field correction");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nx; ii++)
        fprintf(fid, "%.8e %.8e\n", xx[ii], GG[ii]);

    fclose(fid);

    // Output for static density response
    sprintf(out_name, "sdr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
      perror("Error while creating the output file for the static density response");
      exit(EXIT_FAILURE);
    }
    double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
    double ff = 4*lambda*in.rs/M_PI;
    double sdr;
    for (int ii=0 ; ii<in.nx; ii++){
	sdr = -(3.0/2.0)*in.Theta*phi[idx2(ii,0,in.nx)]/
	  (1.0 + ff/(xx[ii]*xx[ii])*(1.0 - GG[ii])*phi[idx2(ii,0,in.nx)]);
	fprintf(fid, "%.8e %.8e\n", xx[ii], sdr);
      }
    fclose(fid);

    // Output for ideal Lindhard density response
    sprintf(out_name, "idr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
      perror("Error while creating the output file for the ideal density response");
      exit(EXIT_FAILURE);
    }
    for (int ii=0; ii<in.nx; ii++){
      for (int jj=0; jj<in.nl; jj++){
        fprintf(fid, "%.8e ", phi[idx2(ii,jj,in.nx)]);
      }
      fprintf(fid,"\n");
    }
    fclose(fid);

    // Output for static structure factor in the Hartree-Fock approximation
    sprintf(out_name, "ssfHF_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static structure factor (HF)");
        exit(EXIT_FAILURE);
    }
    for (int ii = 0; ii < in.nx; ii++)
        fprintf(fid, "%.8e %.8e\n", xx[ii], SSHF[ii]);

    fclose(fid);

    // Output for the interaction energy
    sprintf(out_name, "uint_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
    fid = fopen(out_name, "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the interaction energy");
        exit(EXIT_FAILURE);
    }
    fprintf(fid, "%.8e\n", compute_uex(SS, xx, in));
    fclose(fid);

}


// write binary file to use as initial guess (or restart)
void write_guess(double *SS, double *GG, input in){

  // Name of output file
  char out_name[100];
  sprintf(out_name, "restart_rs%.3f_theta%.3f_%s.bin", in.rs, in.Theta, in.theory);

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(out_name, "wb");
  if (fid == NULL) {
    fprintf(stderr,"Error while creating file for restart");
    exit(EXIT_FAILURE);
  }

  // Static structure factor 
  fwrite(&in, sizeof(input), 1, fid);

  // Static structure factor 
  fwrite(SS, sizeof(double), in.nx, fid);

  // Static local field correction
  fwrite(GG, sizeof(double), in.nx, fid);

  // Close binary file
  fclose(fid);

}


// read binary file to use as initial guess (or restart)
void read_guess(double *SS, double *GG, input in){

  // Variables
  input in_load;

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(in.guess_file, "rb");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file with density response");
    exit(EXIT_FAILURE);
  }

  // Check that the data for the guess file is consistent
  fread(&in_load, sizeof(input), 1, fid);
  if (in_load.nx != in.nx || in_load.dx != in.dx || in_load.xmax != in.xmax){
    fprintf(stderr,"Grid from guess file is incompatible with input\n");
    fclose(fid);
    exit(EXIT_FAILURE);
  }
  
  // Static structure factor in the Hartree-Fock approximation
  fread(SS, sizeof(double), in_load.nx, fid);

  // Static structure factor in the Hartree-Fock approximation
  fread(GG, sizeof(double), in_load.nx, fid);

  // Close binary file
  fclose(fid);
	    
}
