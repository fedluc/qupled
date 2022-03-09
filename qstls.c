#include <omp.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <string.h>
#include "solvers.h"
#include "utils.h"
#include "stls.h"
#include "qstls.h"

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE QSTLS EQUATIONS
// -------------------------------------------------------------------

void solve_qstls(input in, bool verbose) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *SS_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;

  // Arrays for QSTLS solution
  double *psi = NULL;
  double *psi_fixed = NULL;

  // Limit calculations to finite temperatures
  if (in.Theta == 0) {
    fprintf(stderr, "Zero temperature calculations are not "
	    "implemented for the quantum schemes\n");
    exit(EXIT_FAILURE);
  }
  
  // Allocate arrays
  alloc_stls_arrays(in, &xx, &phi, &GG, &SS_new, &SS, &SSHF);
  alloc_qstls_arrays(in, &psi, &psi_fixed);

  // Initialize arrays that are not modified by the iterative procedure
  init_fixed_stls_arrays(&in, xx, phi, SSHF, verbose);
  init_fixed_qstls_arrays(psi_fixed, xx, in, verbose);
			  
  // Initial guess
  initial_guess_qstls(xx, SS, SSHF, psi, phi, in);

  // Iterative procedure
  qstls_iterations(SS, SS_new, SSHF, psi, psi_fixed,
		   phi, xx, in, verbose);
  
  // Internal energy
  if (verbose) printf("Internal energy: %.10f\n",
		      compute_internal_energy(SS, xx, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text_qstls(SS, psi, phi, SSHF, xx, in);
  write_guess_qstls(SS, psi, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_stls_arrays(xx, phi, GG, SS_new, SS, SSHF);
  free_qstls_arrays(psi, psi_fixed);

}


// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_qstls_arrays(input in, double **psi, double **psi_fixed){

  *psi = malloc( sizeof(double) * in.nx * in.nl);
  if (*psi == NULL) {
    fprintf(stderr, "Failed to allocate memory for the auxiliary density response\n");
    exit(EXIT_FAILURE);
  }
  
    
  *psi_fixed = malloc( sizeof(double) * in.nx * in.nl * in.nx);
  if (*psi_fixed == NULL) {
    fprintf(stderr, "Failed to allocate memory for the fixed component of the "
	    "auxiliary density response\n");
    exit(EXIT_FAILURE);
  }

  
}

void free_qstls_arrays(double *psi, double *psi_fixed){

  free(psi);
  free(psi_fixed);
 
}

// ---------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE INITIAL GUESS
// ---------------------------------------------------------------------

void initial_guess_qstls(double *xx, double *SS, double *SSHF,
			 double *psi, double *phi, input in){

  if (strcmp(in.qstls_guess_file,"NO_FILE")==0){

    // Auxilirary density response
    for (int ii=0; ii<in.nx; ii++){
      for (int ll=0; ll<in.nl; ll++){
	psi[idx2(ii,ll,in.nx)] = 0.0;
      }
    }

    // Static structure factor
    compute_ssf_qstls(SS, SSHF, psi, phi, xx, in);
    
  }
  else {

    // Read from file
    read_guess_qstls(SS, psi, in);
    
  }
  
}

// -------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE QSTLS SCHEME
// -------------------------------------------------------------------

void qstls_iterations(double *SS, double *SS_new,
		      double *SSHF, double *psi,
		      double *psi_fixed, double *phi,
		      double *xx, input in,
		      bool verbose){

  double iter_err = 1.0;
  int iter_counter = 0;
  
  if (verbose) printf("SSF calculation...\n");

  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    double tic = omp_get_wtime();
    
    // Update auxiliary density response
    compute_adr(psi, psi_fixed, SS, xx, in);
    
    // Update static structure factor
    compute_ssf_qstls(SS_new, SSHF, psi, phi, xx, in);
    
    // Update diagnostic
    iter_counter++;
    iter_err = stls_err(SS, SS_new, in);
    stls_update(SS, SS_new, in);
      
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

// ------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FIXED COMPONENT OF THE AUXILIARY RESPONSE
// ------------------------------------------------------------------------

// Initialize QSTLS arrays that do not depend on the iterations
void init_fixed_qstls_arrays(double *psi_fixed, double *xx,
			     input in, bool verbose){

  if (strcmp(in.qstls_fixed_file,"NO_FILE")==0){

    // Compute fixed component of the auxiliary density response and
    // store to file
    if (verbose) printf("Fixed component of the QSTLS auxiliary "
			"density response function: ");
    
    compute_adr_fixed(psi_fixed, xx, in);
    write_fixed_qstls(psi_fixed, in);
    
    if (verbose) printf("Done.\n");
    
  }
  else {

    // Read from file
    read_fixed_qstls(psi_fixed, in);
    
  }
  
}


struct adr_fixed_part1_params {

  double mu;
  double Theta;
  gsl_spline *psi_fixed_part1_sp_ptr;
  gsl_interp_accel *psi_fixed_part1_acc_ptr;

};

struct adr_fixed_part2_params {

  double Theta;
  double ww;
  double xx;
  double ll;
  double qq;

};

// Fixed component of the auxiliary density response
void compute_adr_fixed(double *psi_fixed, double *xx, input in) {

  // Parallel calculations
  #pragma omp parallel
  {  

    double err;
    size_t nevals;
    double xx2, xw, tmax, tmin;
    double *psi_fixed_part1  = malloc( sizeof(double) * in.nx);
    if (psi_fixed_part1 == NULL){
      fprintf(stderr, "Failed to allocate memory for fixed component of the "
	      "the auxiliary density response\n");
      exit(EXIT_FAILURE);
    }
      
    // Declare accelerator and spline objects
    gsl_spline *psi_fixed_part1_sp_ptr;
    gsl_interp_accel *psi_fixed_part1_acc_ptr;
    
    // Allocate the accelerator and the spline objects
    psi_fixed_part1_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    psi_fixed_part1_acc_ptr = gsl_interp_accel_alloc();
    
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);

    // Loop over xx (wave-vector)
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){    
      
      // Loop over ll (Matsubara frequencies)
      for (int ll=0; ll<in.nl; ll++){
	
	// Integration function
	gsl_function ff_int_part1, ff_int_part2;
	if (ll == 0){
	  ff_int_part1.function = &adr_fixed_part1_partial_xw0;
	  ff_int_part2.function = &adr_fixed_part2_partial_xwq0;
	}
	else {
	  ff_int_part1.function = &adr_fixed_part1_partial_xwl;
	  ff_int_part2.function = &adr_fixed_part2_partial_xwql;
	}
	
	// Loop over w (wave-vector)
	for (int jj=0; jj<in.nx; jj++){
	  
	  // Integration limits for the integration over t
	  xx2 = xx[ii]*xx[ii];
	  xw = xx[ii]*xx[jj];
	  tmin = xx2 - xw;
	  tmax = xx2 + xw;

	  // Construct integrand for the integral over q
	  for (int kk=0; kk<in.nx; kk++) {
	    
	    if (xx[kk] > 0.0){
	      
	      // Integration over t
	      struct adr_fixed_part2_params ppart2 = {in.Theta,xx[jj],xx[ii],ll,xx[kk]};
	      ff_int_part2.params = &ppart2;
	      gsl_integration_cquad(&ff_int_part2,
				    tmin, tmax,
				    0.0, QUAD_REL_ERR,
				    wsp,
				    &psi_fixed_part1[kk], &err, &nevals);
	    }
	    
	    else psi_fixed_part1[kk] = 0.0;
	    
	  }
	  gsl_spline_init(psi_fixed_part1_sp_ptr, xx, psi_fixed_part1, in.nx);
	  
	  // Integral over q
	  struct adr_fixed_part1_params ppart1 = {in.mu,in.Theta,psi_fixed_part1_sp_ptr,psi_fixed_part1_acc_ptr};
	  ff_int_part1.params = &ppart1;
	  gsl_integration_cquad(&ff_int_part1,
				xx[0], xx[in.nx-1],
				0.0, QUAD_REL_ERR,
				wsp,
				&psi_fixed[idx3(ii, ll, jj, in.nx, in.nl)],
				&err, &nevals);
	}
      }
    }

    // Free memory
    free(psi_fixed_part1);
    gsl_integration_cquad_workspace_free(wsp);
    gsl_spline_free(psi_fixed_part1_sp_ptr);
    gsl_interp_accel_free(psi_fixed_part1_acc_ptr);
    
  }
  
}


// Partial auxiliary density response (vectors = {x,w}, frequency = l)
double adr_fixed_part1_partial_xwl(double qq, void* pp) {
  
  struct adr_fixed_part1_params* params = (struct adr_fixed_part1_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  gsl_spline* psi_fixed_part1_sp_ptr = (params->psi_fixed_part1_sp_ptr);
  gsl_interp_accel* psi_fixed_part1_acc_ptr = (params->psi_fixed_part1_acc_ptr);
  double qq2 = qq*qq;
  double fft = gsl_spline_eval(psi_fixed_part1_sp_ptr, qq, psi_fixed_part1_acc_ptr);

  return qq/(exp(qq2/Theta - mu) + 1.0)*fft;

}


// Partial auxiliary density response (vectors = {x,w}, frequency = 0)
double adr_fixed_part1_partial_xw0(double qq, void* pp) {

  struct adr_fixed_part1_params* params = (struct adr_fixed_part1_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  gsl_spline* psi_fixed_part1_sp_ptr = (params->psi_fixed_part1_sp_ptr);
  gsl_interp_accel* psi_fixed_part1_acc_ptr = (params->psi_fixed_part1_acc_ptr);
  double qq2 = qq*qq;
  double fft = gsl_spline_eval(psi_fixed_part1_sp_ptr, qq, psi_fixed_part1_acc_ptr);

  return qq/(exp(qq2/Theta - mu) + exp(-qq2/Theta + mu) + 2.0)*fft;

}


// Partial auxiliary density response (vectors = {x,q,w}, frequency = l)
double adr_fixed_part2_partial_xwql(double tt, void* pp) {
  
  struct adr_fixed_part2_params* params = (struct adr_fixed_part2_params*)pp;
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

// Partial auxiliary density response (vectors = {x,q,w}, frequency = 0)
double adr_fixed_part2_partial_xwq0(double tt, void* pp) {

  struct adr_fixed_part2_params* params = (struct adr_fixed_part2_params*)pp;
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

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE CHANGING COMPONENT OF THE AUXILIARY RESPONSE
// ---------------------------------------------------------------------------

struct adr_params {

  gsl_spline *adr_fixed_sp_ptr;
  gsl_interp_accel *adr_fixed_acc_ptr;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};


void compute_adr(double *psi, double *psi_fixed, double *SS, 
		 double *xx, input in){
  
  // Parallel calculations
  #pragma omp parallel
  {  

    double err;
    size_t nevals;
    double norm_fact;
    double *adr_fixed  = malloc( sizeof(double) * in.nx);
    if (adr_fixed == NULL){
      fprintf(stderr, "Failed to allocate memory for fixed component of the "
	      "the auxiliary density response\n");
      exit(EXIT_FAILURE);
    }
    
    // Declare accelerator and spline objects
    gsl_spline *ssf_sp_ptr;
    gsl_interp_accel *ssf_acc_ptr;
    gsl_spline *adr_fixed_sp_ptr;
    gsl_interp_accel *adr_fixed_acc_ptr;
  
    // Allocate the accelerator and the spline objects
    adr_fixed_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    adr_fixed_acc_ptr = gsl_interp_accel_alloc();
    ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    ssf_acc_ptr = gsl_interp_accel_alloc();
    
    // Interpolate SSF
    gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);
    
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);
    
    // Integration function
    gsl_function ff_int;
    ff_int.function = &adr_partial;
    
    #pragma omp for // Distribute for loop over the threads
    for (int ii=0; ii<in.nx; ii++){
     for (int ll=0; ll<in.nl; ll++){
  
       // Interpolate solution of q-t integration
       for (int jj=0; jj<in.nx; jj++){
	 adr_fixed[jj] = psi_fixed[idx3(ii,ll,jj,in.nx,in.nl)];
       }
       gsl_spline_init(adr_fixed_sp_ptr, xx, adr_fixed, in.nx);

       // Integral over w 
       struct adr_params pp = {adr_fixed_sp_ptr,adr_fixed_acc_ptr,
				       ssf_sp_ptr, ssf_acc_ptr};
       ff_int.params = &pp;
       gsl_integration_cquad(&ff_int,
			     xx[0], xx[in.nx-1],
			     0.0, QUAD_REL_ERR,
			     wsp,
			     &psi[idx2(ii,ll,in.nx)], 
			     &err, &nevals);
       
       // Assign output
       if (ll == 0) norm_fact = -3.0/(4.0*in.Theta);
       else norm_fact = -3.0/8.0;
       psi[idx2(ii,ll,in.nx)] *= norm_fact;
       
     }
    }  
    
    // Free memory
    free(adr_fixed);
    gsl_integration_cquad_workspace_free(wsp);
    gsl_spline_free(adr_fixed_sp_ptr);
    gsl_interp_accel_free(adr_fixed_acc_ptr);
    gsl_spline_free(ssf_sp_ptr);
    gsl_interp_accel_free(ssf_acc_ptr);
    
  }

}


double adr_partial(double ww, void* pp) {
  
  struct adr_params* params = (struct adr_params*)pp;
  gsl_spline* adr_fixed_sp_ptr = (params->adr_fixed_sp_ptr);
  gsl_interp_accel* adr_fixed_acc_ptr = (params->adr_fixed_acc_ptr);
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);

  return ww*gsl_spline_eval(adr_fixed_sp_ptr, ww, adr_fixed_acc_ptr)
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
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

// write text files for output
void write_text_qstls(double *SS, double *psi, double *phi, 
		      double *SSHF, double *xx, input in){

  bool finite_temperature = false;

  // Check the value of the temperature
  if (in.Theta > 0) finite_temperature = true;
  
  // Static structure factor
  write_text_ssf(SS, xx, in);

  // Static structure factor within the Hartree-Fock approximation
  write_text_ssf_HF(SSHF, xx, in);
  
  // Static local field correction
  if (finite_temperature)
    write_text_slfc_qstls(psi, phi, xx, in);
  
  // Static density response
  if (finite_temperature)
    write_text_sdr_qstls(psi, phi, xx, in);
  
  // Ideal density response
  if (finite_temperature)
    write_text_idr(phi, in);

  // Auxiliary density response
  if (finite_temperature)
    write_text_adr(psi, in);

  // Radial distribution function
  write_text_rdf(SS, xx, in);

  // Interaction energy
  write_text_uint(SS, xx, in);

}

// write static local field correction to text file
void write_text_slfc_qstls(double *psi, double *phi, double *xx,  input in){

  FILE* fid;
  char out_name[100];
  
  sprintf(out_name, "slfc_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the static local field correction\n");
    exit(EXIT_FAILURE);
  }

  for (int ii = 0; ii < in.nx; ii++)
    fprintf(fid, "%.8e %.8e\n", xx[ii], psi[idx2(ii,0,in.nx)]/phi[idx2(ii,0,in.nx)]);
  
  fclose(fid);
  
}


// write static density response to text file
void write_text_sdr_qstls(double *psi, double *phi, double *xx,  input in){

  FILE* fid;
  char out_name[100];
  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double ff = 4*lambda*in.rs/M_PI;
  double sdr;
  
  sprintf(out_name, "sdr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the static density response");
    exit(EXIT_FAILURE);
  }
  
  for (int ii=0 ; ii<in.nx; ii++){
    sdr = -(3.0/2.0)*in.Theta*phi[idx2(ii,0,in.nx)]/
      (1.0 + ff/(xx[ii]*xx[ii])*(phi[idx2(ii,0,in.nx)]
				 - psi[idx2(ii,0,in.nx)]));
    fprintf(fid, "%.8e %.8e\n", xx[ii], sdr);
  }

  fclose(fid);
  
}


// write auxiliary density response to text file
void write_text_adr(double *psi, input in){
  
  FILE* fid;
  char out_name[100];
  
  sprintf(out_name, "adr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    fprintf(stderr, "Error while creating the output file for the auxiliary density response");
    exit(EXIT_FAILURE);
  }
  
  for (int ii=0; ii<in.nx; ii++){
    for (int jj=0; jj<in.nl; jj++){
      fprintf(fid, "%.8e ", psi[idx2(ii,jj,in.nx)]);
    }
    fprintf(fid,"\n");
  }

  fclose(fid);
  
}


// write binary file to use as initial guess (or restart)
void write_guess_qstls(double *SS, double *psi, input in){

  // Name of output file
  char out_name[100];
  sprintf(out_name, "restart_rs%.3f_theta%.3f_%s.bin", in.rs, in.Theta, in.theory);

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(out_name, "wb");
  if (fid == NULL) {
    fprintf(stderr,"Error while creating file for initial guess or restart\n");
    exit(EXIT_FAILURE);
  }

  // Input data
  fwrite(&in.nx, sizeof(int), 1, fid);
  fwrite(&in.nl, sizeof(int), 1, fid);
  fwrite(&in.dx, sizeof(double), 1, fid);
  fwrite(&in.xmax, sizeof(double), 1, fid);

  // Static structure factor 
  fwrite(SS, sizeof(double), in.nx, fid);

  // Auxilliary density response 
  fwrite(psi, sizeof(double), in.nx*in.nl, fid);

  // Close binary file
  fclose(fid);

}


// read binary file to use as initial guess (or restart)
void read_guess_qstls(double *SS, double *psi, input in){

  // Variables
  size_t it_read;
  int nx_file;
  int nl_file;
  double dx_file;
  double xmax_file;

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(in.qstls_guess_file, "rb");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file for initial guess or restart\n");
    exit(EXIT_FAILURE);
  }

  // Initialize number of items read from input file
  it_read = 0;

  // Check that the data for the guess file is consistent
  it_read += fread(&nx_file, sizeof(int), 1, fid);
  it_read += fread(&nl_file, sizeof(int), 1, fid);
  it_read += fread(&dx_file, sizeof(double), 1, fid);
  it_read += fread(&xmax_file, sizeof(double), 1, fid);
  check_guess_qstls(nx_file, dx_file, xmax_file, nl_file, in.Theta,
		   in, it_read, 4, fid, true, true, false);
  
  
  // Static structure factor
  it_read += fread(SS, sizeof(double), nx_file, fid);
  
  // Auxilliary density response
  it_read += fread(psi, sizeof(double), nx_file * nl_file, fid);

  // Check that all the expected items where read
  check_guess_qstls(nx_file, dx_file, xmax_file, nl_file, in.Theta,
		   in, it_read, nx_file + nl_file*nx_file + 4,
		   fid, false, true, true);
  

  // Close binary file
  fclose(fid);
	    
}

// write binary file to store the fixed component of the auxilliary density response
void write_fixed_qstls(double *psi_fixed, input in){

  // Name of output file
  char out_name[100];
  sprintf(out_name, "fixed_rs%.3f_theta%.3f_QSTLS.bin", in.rs, in.Theta);

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(out_name, "wb");
  if (fid == NULL) {
    fprintf(stderr,"Error while creating file for the fixed component of the auxilliary density response\n");
    exit(EXIT_FAILURE);
  }

  // Input data
  fwrite(&in.nx, sizeof(int), 1, fid);
  fwrite(&in.nl, sizeof(int), 1, fid);
  fwrite(&in.dx, sizeof(double), 1, fid);
  fwrite(&in.xmax, sizeof(double), 1, fid);
  fwrite(&in.Theta, sizeof(double), 1, fid);

  // Fixed component of the auxiliary density response
  fwrite(psi_fixed, sizeof(double), in.nx * in.nl * in.nx, fid);

  // Close binary file
  fclose(fid);

}

// read binary file with the fixed component of the auxilliary density response
void read_fixed_qstls(double *psi_fixed, input in){

  // Variables
  size_t it_read;
  int nx_file;
  int nl_file;
  double dx_file;
  double xmax_file;
  double Theta_file;

  // Open binary file
  FILE *fid = NULL;
  fid = fopen(in.qstls_fixed_file, "rb");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file for initial guess or restart\n");
    exit(EXIT_FAILURE);
  }

  // Initialize number of items read from input file
  it_read = 0;

  // Check that the data for the guess file is consistent
  it_read += fread(&nx_file, sizeof(int), 1, fid);
  it_read += fread(&nl_file, sizeof(int), 1, fid);
  it_read += fread(&dx_file, sizeof(double), 1, fid);
  it_read += fread(&xmax_file, sizeof(double), 1, fid);
  it_read += fread(&Theta_file, sizeof(double), 1, fid);
  check_guess_qstls(nx_file, dx_file, xmax_file, nl_file, Theta_file,
		   in, it_read, 5, fid, true, true, false);
  

  // Fixed component of the auxiliary density response
  it_read += fread(psi_fixed, sizeof(double), nx_file * nl_file * nx_file, fid);

  // Check that all items where read and the end-of-file was reached
  check_guess_qstls(nx_file, dx_file, xmax_file, nl_file, Theta_file,
		   in, it_read, nx_file*nl_file*nx_file + 5, fid,
		   false, true, true);
  
  // Close binary file
  fclose(fid);
	    
}

// Check consistency of the guess data
void check_guess_qstls(int nx, double dx, double xmax, int nl,
		       double Theta, input in, size_t it_read,
		       size_t it_expected, FILE *fid, bool check_grid,
		       bool check_items, bool check_eof){
  
  int buffer;
  
  // Check that the grid in the guess data is consistent with input
  if (check_grid) {
    
    if (nx != in.nx || fabs(dx-in.dx) > DBL_TOL || fabs(xmax-in.xmax) > DBL_TOL){
      fprintf(stderr,"Grid from guess file is incompatible with input\n");
      fprintf(stderr,"Grid points (nx) : %d (input), %d (file)\n", in.nx, nx);
      fprintf(stderr,"Resolution (dx)  : %.16f (input), %.16f (file)\n", in.dx, dx);
      fprintf(stderr,"Cutoff (xmax)    : %.16f (input), %.16f (file)\n", in.xmax, xmax);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
    if (nl != in.nl){
      fprintf(stderr,"Number of Matsubara frequencies from fixed solution file is incompatible with input\n");
      fprintf(stderr,"Matsubara frequencies (nl) : %d (input), %d (file)\n", in.nl, nl);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
    if (fabs(Theta-in.Theta) > DBL_TOL){
      fprintf(stderr,"Quantum degeneracy parameter from fixed solution file is incompatible with input\n");
      fprintf(stderr,"Degeneracy parameter (theta) : %f (input), %f (file)\n", in.Theta, Theta);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
    
  }

  // Check that all the expected items where read
  if (check_items) {
    if (it_read != it_expected ) {
      fprintf(stderr,"Error while reading file for initial guess or restart.\n");
      fprintf(stderr,"%ld Elements expected, %ld elements read\n", it_read, it_expected);
      fclose(fid);
      exit(EXIT_FAILURE);
    }
  }
  
  // Check for end of file
  if (check_eof){
    it_read = fread(&buffer, sizeof(int), 1, fid); // Trigger end-of-file activation
    if (!feof(fid)) {
      fprintf(stderr,"Error while reading file for initial guess or restart.\n");
      fprintf(stderr,"Expected end of file, but there is still data left to read.\n");
      fclose(fid);
      exit(EXIT_FAILURE);
    }
  }
  
}
