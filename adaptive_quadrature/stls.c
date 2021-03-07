#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "stls.h"


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS EQUATIONS
// -------------------------------------------------------------------

void solve_stls(input in, bool verbose, 
		double **xx_out, double **SS_out, 
		double **SSHF_out, double **GG_out, 
		double **GG_new_out, double **phi_out) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *GG_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;

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
  double iter_err = 1.0;
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
      printf("Elapsed time: %f seconds\n", ((double)toc - (double)tic) / CLOCKS_PER_SEC);
      printf("Residual error: %.5e\n", iter_err);
    }
  }
  if (verbose) printf("Done.\n");
  
  // Internal energy
  if (verbose) printf("Internal energy: %f\n",compute_internal_energy(SS, xx, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text(SS, GG, xx, in);
  if (init_flag) write_bin(phi, SSHF, in); 
  if (verbose) printf("Done.\n");

  // Output to variable or free memory
  if (xx_out != NULL) {
    *xx_out = xx;
    *SS_out = SS;
    *SSHF_out = SSHF;
    *GG_out = GG;
    *GG_new_out = GG_new;
    *phi_out = phi;
  }
  else{
    free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
  }
 
 
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

void free_stls_arrays(double *xx, double *phi, 
		      double *GG, double *GG_new, 
		      double *SS, double *SSHF){

  free(xx);
  free(phi);
  free(SSHF);
  free(SS);
  free(GG);
 
}


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

struct nc_params {

  double Theta;

};

double compute_mu(input in) {
  
  // Variables
  int max_iter = 100;
  double mu_lo = -10.0;
  double mu_hi = 10.0;
  double mu;
  int status, iter;

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
  return mu;

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

void wave_vector_grid(double *xx, input in){
  
  //xx[0] = in.dx;
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
  double plt = M_PI*ll*in.Theta, plt2 = plt*plt;

  // Integration workspace 
  gsl_integration_cquad_workspace *wsp 
    = gsl_integration_cquad_workspace_alloc(1000);

  // Integration function
  gsl_function ff_int;
  if (ll == 0) ff_int.function = &phix0;
  else ff_int.function = &phixl;

  // Normalized ideal Lindhard density 
  for (int ii = 0; ii < in.nx; ii++) {
    
    if (xx[ii] >= 0.01){
      struct phixl_params phixlp = {xx[ii], in.mu, in.Theta, ll};
      ff_int.params = &phixlp;
      gsl_integration_cquad(&ff_int, 
			    xx[0], xx[in.nx-1], 
			    0.0, in.err_min_int, 
			    wsp, 
			    &phil[ii], &err, &nevals);
    }
    else {
      
      if (ll ==0) phil[ii] = sqrt(in.Theta)*gsl_sf_gamma(0.5)
		    *gsl_sf_fermi_dirac_mhalf(in.mu)/2.0;

      else phil[ii] = xx[ii]*xx[ii]/(3*plt2);

    }
  }

  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  
  
}

double phixl(double yy, void* pp) {

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

double phix0(double yy, void* pp) {

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

struct ssfHF_params {

  double xx;
  double mu;
  double Theta;

};

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

void compute_ssfHF(double *SS,  double *xx,  input in){

  double err;
  size_t nevals;

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(in.nx);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &ssfHF;

  // Static structure factor in the Hartree-Fock approximation
  for (int ii = 0; ii < in.nx; ii++) {


    struct ssfHF_params ssfHFp = {xx[ii], in.mu, in.Theta};
    ff_int.params = &ssfHFp;
    gsl_integration_cquad(&ff_int,
			  xx[0], xx[in.nx-1],
			  0.0, in.err_min_int,
			  wsp,
			  &SS[ii], &err, &nevals);

    SS[ii] += 1.0;

  }
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);

}

double csch2(double x){
  double csch = 1.0/sinh(x);
  return csch*csch;
}

double coth(double x){
  return 1.0/tanh(x);
}


void compute_AA(double *AA, double *xx,  input in){

  double fact1 = 2.0/(9.0*in.Theta*in.Theta);
  double xx2, yy;
  
  for (int ii=0; ii<in.nx; ii++){
    xx2 = xx[ii]*xx[ii];
    yy = xx2/(2.0*in.Theta);
    AA[ii] = fact1*(csch2(yy) + coth(yy)/yy);
  }
}

double Axl2(double xx, int ll, input in){
  
  double xx2 = xx*xx;
  double xx4 = xx2*xx2;
  double tplT = 2*M_PI*ll*in.Theta;
  double tplT2 = tplT*tplT;
  double Axl = (4.0/3.0)*xx2/(tplT2 + xx4);
  
  return Axl*Axl;

}


void compute_ssf(double *SS, double *SSHF, double *GG, 
		 double *phi, double *xx, input in){

  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  /* double ff = 4.0*lambda*in.rs/M_PI; */
  /* double ff3_2T = 3.0*in.Theta*ff/2.0; */
  double pilambda = M_PI*lambda;
  double ff = 4*lambda*lambda*in.rs;
  double ff3_2T = 3.0*in.Theta*ff/2.0;
  double xx2, BB, BB_tmp, phixl;
  for (int ii=0; ii<in.nx; ii++){

    xx2 = xx[ii]*xx[ii];
    BB = 0.0;
    if ( xx[ii] > 0.0){

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
    = gsl_integration_cquad_workspace_alloc(in.nx);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &slfc;

  // Static local field correction
  for (int ii = 0; ii < in.nx; ii++) {
    struct slfc_params slfcp = {xx[ii], ssf_sp_ptr, ssf_acc_ptr};
    ff_int.params = &slfcp;
    gsl_integration_cquad(&ff_int,
  			  xx[0], xx[in.nx-1],
  			  0.0, in.err_min_int,
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
        return yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0)
          * (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(xx - yy)));
      }
      else if (xx < yy) {
        return yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0)
          * (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(yy - xx)));
      }
      else {
        return yy2 * (gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1.0);
      }

    }
    else {
      return 0;
    }



}



// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

struct uex_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};


double compute_internal_energy(double *SS, double *xx,  input in) {

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
    = gsl_integration_cquad_workspace_alloc(in.nx);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &uex;

  // Internal energy
  struct uex_params uexp = {ssf_sp_ptr, ssf_acc_ptr};
  ff_int.params = &uexp;  
  gsl_integration_cquad(&ff_int,
			xx[0], xx[in.nx-1],
			0.0, in.err_min_int,
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


// write text files with SSF and SLFC
void write_text(double *SS, double *GG, double *xx, input in){


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
void write_bin(double *phi, double *SSHF, input in){


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
  fwrite(phi, sizeof(double), in.nx*in.nl, fid);

  // Static structure factor in the Hartree-Fock approximation
  fwrite(SSHF, sizeof(double), in.nx, fid);

  // Close binary file
  fclose(fid);

}

// read text file with SSF and SLFC
void read_text(double *SS, double *GG, double *xx, input in){
  
  /* FILE *fid; */

  /* // Read files once to get the size  */
  /* fid = fopen(filename, "r"); */

  /* if (fid == NULL) { */
  /*   perror("Error while opening ssf file"); */
  /*   exit(EXIT_FAILURE); */
  /* } */
 
  
 

}
 

// read binary file with density response information
void read_bin(input *in, double **xx, double **phi, 
	      double **GG, double **GG_new, 
	      double **SS, double **SSHF){

  // Variables
  double *xx_local = NULL; 
  double *phi_local = NULL;
  double *GG_local = NULL;
  double *GG_new_local = NULL;
  double *SS_local = NULL;
  double *SSHF_local = NULL;
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
  fread(phi_local, sizeof(double), in_load.nx*in_load.nl, fid);

  // Static structure factor in the Hartree-Fock approximation
  fread(SSHF_local, sizeof(double), in_load.nx, fid);

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
