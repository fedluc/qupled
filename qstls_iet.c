#include <omp.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <string.h>
#include "solvers.h"
#include "stls.h"
#include "qstls.h"
#include "qstls_iet.h"
#include "stls_iet.h"

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE QSTLS EQUATIONS
// -------------------------------------------------------------------

void solve_qstls_iet(input in, bool verbose) {

  // Arrays for STLS solution
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *SS_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;

  // Arrays for QSTLS solution
  double *psi = NULL;
  double *psi_new = NULL;
  double *psi_fixed_qstls = NULL;
  double *bf = NULL;

  // Allocate arrays
  // Note: GG is not needed for QSTLS, but we keep it here so that
  // we can reuse some stls routines 
  alloc_stls_arrays(in, &xx, &phi, &GG, &SS_new, &SS, &SSHF);
  psi = malloc( sizeof(double) * in.nx * in.nl);  
  psi_new = malloc( sizeof(double) * in.nx * in.nl);  
  psi_fixed_qstls = malloc( sizeof(double) * in.nx * in.nl * in.nx);  
  bf = malloc(sizeof(double) * in.nx);

  // Initialize STLS arrays that are not modified by the iterative procedure
  init_fixed_stls_arrays(&in, xx, phi, SSHF, verbose);
  compute_bridge_function(bf, xx, in);

  // Initial guess
  if (strcmp(in.qstls_guess_file,"NO_FILE")==0){
    for (int ii=0; ii<in.nx; ii++){
      for (int ll=0; ll<in.nl; ll++){
  	psi[idx2(ii,ll,in.nx)] = 0.0;
      }
    }
    compute_ssf_qstls_iet(SS, SSHF, psi, phi, bf, xx, in);
  }
  else {
    read_guess_qstls(SS, psi, in);
  }

  // Initialize QSTLS arrays that are not modified by the iterative procedure
  if (verbose) printf("Fixed component of the qstls auxiliary response function: ");
  if (strcmp(in.qstls_fixed_file,"NO_FILE")==0){
    compute_adr_fixed(psi_fixed_qstls, xx, in);
    write_fixed_qstls(psi_fixed_qstls, in);
  }
  else {
    read_fixed_qstls(psi_fixed_qstls, in);
  }
  if (verbose) printf("Done.\n");

  // Initialize QSTLS-IET arrays that are not modified by the iterative procedure
  if (verbose) printf("Fixed component of the auxiliary response function: ");
  if (strcmp(in.qstls_iet_fixed_file,"NO_FILE")==0){
    compute_adr_iet_fixed(xx, in);
  }
  if (verbose) printf("Done.\n");
 
  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF calculation...\n");
  double iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    double tic = omp_get_wtime();
    
    // Update auxiliary function
    compute_adr_iet(psi_new, psi, psi_fixed_qstls, phi, SS, bf, xx, in);

    // Update SSF
    compute_ssf_qstls_iet(SS_new, SSHF, psi_new, phi, bf, xx, in);

    // Prepare auxiliary function for next iteration
    for (int ii=0; ii<in.nx; ii++){
      for (int ll=0; ll<in.nl; ll++){
    	psi[idx2(ii,ll,in.nx)] = psi_new[idx2(ii,ll,in.nx)];
      }
    }

    // Update diagnostic
    iter_err = 0.0;
    iter_counter++;
    for (int ii=0; ii<in.nx; ii++) {
      iter_err += (SS_new[ii] - SS[ii]) * (SS_new[ii] - SS[ii]);
      SS[ii] = in.a_mix*SS_new[ii] + (1-in.a_mix)*SS[ii];
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
  if (verbose) printf("Internal energy: %.10f\n",compute_internal_energy(SS, xx, in));
  
  // Output to file
  if (verbose) printf("Writing output files...\n");
  write_text_qstls_iet(SS, psi, phi, SSHF, bf, xx, in);
  write_guess_qstls(SS, psi, in);
  write_bridge_function(bf, xx, in);
  if (verbose) printf("Done.\n");

  // Free memory
  free_stls_arrays(xx, phi, GG, SS_new, SS, SSHF);
  free(psi);
  free(psi_new);
  free(psi_fixed_qstls);
  free(bf);

}


// ------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FIXED COMPONENT OF THE AUXILIARY RESPONSE
// ------------------------------------------------------------------------

struct adr_iet_fixed_params {

  double mu;
  double Theta;
  double ll;
  double xx;
  double uu;
  double ww;

};

// Fixed component of the auxiliary density response within the QSTLS-IET scheme
void compute_adr_iet_fixed(double *xx, input in) {

  // Parallel calculations
  #pragma omp parallel
  {  

    double err;
    size_t nevals;
    double adr_iet_fixed;
     
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);

    // Loop over xx (wave-vector)
    #pragma omp for // Distribute loops over the threads
    for (int ii=0; ii<in.nx; ii++){

      // Open binary file for output
      char out_name[100];
      sprintf(out_name, "psi_fixed_theta%.3f_xx%.5f.bin", in.Theta, xx[ii]);
      FILE *fid = NULL;
      fid = fopen(out_name, "wb");
      if (fid == NULL) {
	fprintf(stderr,"Error while creating file for fixed component of the auxilliary response function");
	exit(EXIT_FAILURE);
      }

      // Write input data to file
      fwrite(&in.nx, sizeof(int), 1, fid);
      fwrite(&in.nl, sizeof(int), 1, fid);
      fwrite(&in.dx, sizeof(double), 1, fid);
      fwrite(&in.xmax, sizeof(double), 1, fid);
      fwrite(&in.Theta, sizeof(double), 1, fid);
      
      // Loop over ll (Matsubara frequencies)
      for (int ll=0; ll<in.nl; ll++){
	
    	// Integration function
    	gsl_function ff_int;
    	if (ll == 0){
    	  ff_int.function = &adr_iet_fixed_partial_xuw0;
    	}
    	else {
    	  ff_int.function = &adr_iet_fixed_partial_xuwl;
    	}
	
    	// Loop over u (wave-vector)
    	for (int jj=0; jj<in.nx; jj++){

  	  // Loop over w
  	  for (int kk=0; kk<in.nx; kk++) {
	    
  	    // Integral over y
  	    if (xx[ii] == 0.0 || xx[jj] == 0.0 || xx[kk] == 0.0){
  	      // No need to compute the integral over y in this cases
  	      adr_iet_fixed = 0.0;
  	    }
  	    else{
  	      // Compute the integral and store it
  	      struct adr_iet_fixed_params pp = {in.mu,in.Theta,ll,xx[ii],xx[jj],xx[kk]};
  	      ff_int.params = &pp;
  	      gsl_integration_cquad(&ff_int,
  				    xx[0], xx[in.nx-1],
  				    0.0, QUAD_REL_ERR,
  				    wsp,
  				    &adr_iet_fixed,
  				    &err, &nevals);
           }
	    
	    // Write result to output file
	    fwrite(&adr_iet_fixed, sizeof(double), 1, fid);
	    
	  }
    	}
      }

      // Close binary file for output
      fclose(fid);

    }

    // Free memory
    gsl_integration_cquad_workspace_free(wsp);
 
  }
  
}

// Partial auxiliary density response (vectors = {x,u,w}, frequency = 0)
double adr_iet_fixed_partial_xuwl(double yy, void* pp) {
  
  struct adr_iet_fixed_params* params = (struct adr_iet_fixed_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  double ll = (params->ll);
  double xx = (params->xx);
  double uu = (params->uu);
  double ww = (params->ww);
  double xx2 = xx*xx, ww2 =  ww*ww, uu2 = uu*uu, fxy = 4.0*xx*yy, 
    yy2 = yy*yy, fplT = 4.0*M_PI*ll*Theta, fplT2 = fplT*fplT,
    umwpx = uu2 - ww2 + xx2,
    logarg = ((umwpx + fxy)*(umwpx + fxy) + fplT2)/
             ((umwpx - fxy)*(umwpx - fxy) + fplT2);

  if (xx == 0 || uu == 0 || ww == 0)
    return 0;
  else
    return yy/(exp(yy2/Theta - mu) + 1.0)*log(logarg);

}


// Partial auxiliary density response (vectors = {x,u,w}, frequency = 0)
double adr_iet_fixed_partial_xuw0(double yy, void* pp) {

  struct adr_iet_fixed_params* params = (struct adr_iet_fixed_params*)pp;
  double mu = (params->mu);
  double Theta = (params->Theta);
  double xx = (params->xx);
  double uu = (params->uu);
  double ww = (params->ww);
  double xx2 = xx*xx, ww2 =  ww*ww, uu2 = uu*uu, fxy = 4.0*xx*yy, 
    yy2 = yy*yy, umwpx = uu2 - ww2 + xx2,
    logarg = (umwpx + fxy)/(umwpx - fxy);

  if (logarg < 0.0) logarg = -logarg;

  if (xx == 0 || uu == 0 || ww == 0)
    return 0;
  else
    return yy/(exp(yy2/Theta - mu) + exp(-yy2/Theta + mu) + 2.0)*
      ((yy2 - umwpx*umwpx/(16.0*xx2))*log(logarg) + (yy/xx)*umwpx/2.0);

}

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE CHANGING COMPONENT OF THE AUXILIARY RESPONSE
// ---------------------------------------------------------------------------

struct adr_iet_part1_params {

  gsl_spline *adr_part1_sp_ptr;
  gsl_interp_accel *adr_part1_acc_ptr;

};


struct adr_iet_part2_params {

  gsl_spline *adr_part2_sp_ptr;
  gsl_interp_accel *adr_part2_acc_ptr;

};

void compute_adr_iet(double *psi_new, double *psi, double *psi_fixed_qstls, 
		     double *phi, double *SS, double *bf, double *xx, 
		     input in){


  // QSTLS component of the auxilliary response function
  compute_adr(psi_new, psi_fixed_qstls, SS, xx, in);

  // Parallel calculations
  #pragma omp parallel
  {  

    double err;
    size_t nevals;
    double psi_tmp;
    double wmax, wmin;
    double *adr_part2  = malloc( sizeof(double) * in.nx);
    double *adr_part1  = malloc( sizeof(double) * in.nx);
    
    // Declare accelerator and spline objects
    gsl_spline *adr_part1_sp_ptr;
    gsl_interp_accel *adr_part1_acc_ptr;
    gsl_spline *adr_part2_sp_ptr;
    gsl_interp_accel *adr_part2_acc_ptr;
    
    // Allocate the accelerator and the spline objects
    adr_part1_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    adr_part1_acc_ptr = gsl_interp_accel_alloc();
    adr_part2_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
    adr_part2_acc_ptr = gsl_interp_accel_alloc();
    
    // Integration workspace
    gsl_integration_cquad_workspace *wsp
      = gsl_integration_cquad_workspace_alloc(100);
    
    // Integration function
    gsl_function ff_part1_int, ff_part2_int;
    ff_part1_int.function = &adr_iet_part1_partial;
    ff_part2_int.function = &adr_iet_part2_partial;
    
    // Loop over xx (wave-vector)
    #pragma omp for // Distribute loops over the threads
    for (int ii=0; ii<in.nx; ii++){
      
      // Open binary file with the fixed component of the auxilliary response function
      char out_name[100000];
      size_t it_read = 0;
      int nx_file;
      int nl_file;
      double dx_file;
      double xmax_file;
      double Theta_file;
      FILE *fid = NULL;
      if (strcmp(in.qstls_iet_fixed_file,"NO_FILE")==0){
	sprintf(out_name, "psi_fixed_theta%.3f_xx%.5f.bin", in.Theta, xx[ii]);
      }
      else{
	sprintf(out_name, "%s/psi_fixed_theta%.3f_xx%.5f.bin", in.qstls_iet_fixed_file, in.Theta, xx[ii]);
      }
      fid = fopen(out_name, "rb");
      if (fid == NULL) {
	fprintf(stderr,"Error while reading file for fixed component of the auxilliary response function\n");
	exit(EXIT_FAILURE);
      }

      // Check that the data for the binary file is consistent with input
      it_read += fread(&nx_file, sizeof(int), 1, fid);
      it_read += fread(&nl_file, sizeof(int), 1, fid);
      it_read += fread(&dx_file, sizeof(double), 1, fid);
      it_read += fread(&xmax_file, sizeof(double), 1, fid);
      it_read += fread(&Theta_file, sizeof(double), 1, fid);
      check_guess_qstls(nx_file, dx_file, xmax_file, nl_file, Theta_file,
			in, it_read, 5, fid, true, true, false);
      
      // Loop over the Matsubara frequencies
      for (int ll=0; ll<in.nl; ll++){
	
      	// Loop over u (wave-vector)
      	for (int jj=0; jj<in.nx; jj++){
	  
      	  // Construct integrand over w
      	  for (int kk=0; kk<in.nx; kk++) {
      	    fread(&adr_part2[kk], sizeof(double), 1, fid);
      	    adr_part2[kk] *= xx[kk]*(SS[kk]-1.0);
      	  }
      	  gsl_spline_init(adr_part2_sp_ptr, xx, adr_part2, in.nx);
	  
      	  // Compute integral over w
	  if (xx[ii] == 0.0 || xx[jj] == 0.0){
	    adr_part1[jj] = 0.0;
	  }
	  else {

	    struct adr_iet_part2_params pp_part2 =  {adr_part2_sp_ptr, adr_part2_acc_ptr};
	    wmin = xx[jj] - xx[ii];
	    if (wmin < 0.0) wmin = -wmin;
	    // NOTE: The upper cutoff is at qm - dq for numerical reasons;
	    // The quadrature formula attemps a tiny extrapolation which causes
	    // the interpolation routine to crash.
	    wmax = GSL_MIN(xx[in.nx-2], xx[ii]+xx[jj]);
	    ff_part2_int.params = &pp_part2;
	    gsl_integration_cquad(&ff_part2_int,
				  wmin, wmax,
				  0.0, QUAD_REL_ERR,
				  wsp,
				  &adr_part1[jj], &err, &nevals);

	    // Construct integrand over u ( the -1 is added because the qSTLS contribution is calculated separately)
	    if (in.qstls_iet_static){ 
	      adr_part1[jj] *= (1.0/xx[jj])
		*( (-bf[jj]+1)*SS[jj] - 1 - 
		   (psi[idx2(jj,0,in.nx)]/phi[idx2(jj,0,in.nx)])*(SS[jj]-1));
	    }
	    else {
	      adr_part1[jj] *= (1.0/xx[jj])
		*( (-bf[jj]+1)*SS[jj] - 1 - 
		   (psi[idx2(jj,ll,in.nx)]/phi[idx2(jj,ll,in.nx)])*(SS[jj]-1));
	    }
	    
	  }
	  
      	}
	
      	// Interpolate integrand over u
      	gsl_spline_init(adr_part1_sp_ptr, xx, adr_part1, in.nx);
	
      	// Compute integral over u
      	struct adr_iet_part1_params pp_part1 = {adr_part1_sp_ptr, adr_part1_acc_ptr};
      	ff_part1_int.params = &pp_part1;
      	gsl_integration_cquad(&ff_part1_int,
      			      xx[0], xx[in.nx-1],
      			      0.0, QUAD_REL_ERR,
      			      wsp,
      			      &psi_tmp, &err, &nevals);
	
      	if (ll == 0)
      	  psi_new[idx2(ii,ll,in.nx)] += -3.0/(4.0*in.Theta)*psi_tmp;
      	else
      	  psi_new[idx2(ii,ll,in.nx)] += -(3.0/8.0)*psi_tmp;
	
      }

      // Check that all the binary file was read
      check_guess_qstls(nx_file, dx_file, xmax_file, nl_file, Theta_file,
			in, it_read, 0, fid, false, false, true);
      
      // Close binary file for output
      fclose(fid);
      
    }
    
    // Free memory
    free(adr_part2);
    free(adr_part1);
    gsl_integration_cquad_workspace_free(wsp);
    gsl_spline_free(adr_part1_sp_ptr);
    gsl_interp_accel_free(adr_part1_acc_ptr);
    gsl_spline_free(adr_part2_sp_ptr);
    gsl_interp_accel_free(adr_part2_acc_ptr);
    
  }
  
}


double adr_iet_part1_partial(double uu, void* pp) {

  struct adr_iet_part1_params* params = (struct adr_iet_part1_params*)pp;
  gsl_spline* adr_part1_sp_ptr = (params->adr_part1_sp_ptr);
  gsl_interp_accel* adr_part1_acc_ptr = (params->adr_part1_acc_ptr);
 
  return gsl_spline_eval(adr_part1_sp_ptr, uu, adr_part1_acc_ptr);

}

double adr_iet_part2_partial(double ww, void* pp) {

  struct adr_iet_part2_params* params = (struct adr_iet_part2_params*)pp;
  gsl_spline* adr_part2_sp_ptr = (params->adr_part2_sp_ptr);
  gsl_interp_accel* adr_part2_acc_ptr = (params->adr_part2_acc_ptr);
 
  return gsl_spline_eval(adr_part2_sp_ptr, ww, adr_part2_acc_ptr);
}


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_qstls_iet(double *SS, double *SSHF, double *psi,
			   double *phi, double *bf, double *xx, input in){

  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  double ff = 4*lambda*in.rs/M_PI;
  double xx2, BB, BB_tmp, BB_den, psixl, phixl, bfphixl;

  for (int ii=0; ii<in.nx; ii++){

    if (xx[ii] > 0.0){

      xx2 = xx[ii]*xx[ii];
      BB = 0.0;

      for (int ll=0; ll<in.nl; ll++){

        psixl = psi[idx2(ii,ll,in.nx)];
        phixl = phi[idx2(ii,ll,in.nx)];
	bfphixl = (1.0 - bf[ii])*phixl;
        BB_den = 1.0 + ff/xx2*(1.0 - psixl/bfphixl)*bfphixl;
        BB_tmp = phixl*phixl*(1.0 - psixl/bfphixl)/BB_den;
        if (ll>0) BB_tmp *= 2.0;
        BB += BB_tmp;

      }

      SS[ii] = SSHF[ii] - 3.0/2.0*ff/xx2*in.Theta*(1.0-bf[ii])*BB;

    }
    else
      SS[ii] = 0.0;

  }

}

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------


// write text files for output
void write_text_qstls_iet(double *SS, double *psi, double *phi,
			  double *SSHF, double *bf, double *xx, 
			  input in){


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
    fprintf(fid, "%.8e %.8e\n", xx[ii], psi[idx2(ii,0,in.nx)]/phi[idx2(ii,0,in.nx)]);

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
      (1.0 + ff/(xx[ii]*xx[ii])*((1.0 - bf[ii])*phi[idx2(ii,0,in.nx)] - psi[idx2(ii,0,in.nx)]));
    fprintf(fid, "%.8e %.8e\n", xx[ii], sdr);
  }
  fclose(fid);

  // Output for the auxiliary density response
  sprintf(out_name, "adr_rs%.3f_theta%.3f_%s.dat", in.rs, in.Theta, in.theory);
  fid = fopen(out_name, "w");
  if (fid == NULL) {
    perror("Error while creating the output file for the auxiliary density response");
    exit(EXIT_FAILURE);
  }
  for (int ii=0; ii<in.nx; ii++){
    for (int jj=0; jj<in.nl; jj++){
      fprintf(fid, "%.8e ", psi[idx2(ii,jj,in.nx)]);
    }
    fprintf(fid,"\n");
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
  fprintf(fid, "%.8e\n", compute_internal_energy(SS, xx, in));
  fclose(fid);

}
