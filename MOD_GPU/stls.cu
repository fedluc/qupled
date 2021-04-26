extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "solvers.h"
#include "chemical_potential.h"
#include "stls.h"
}

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS EQUATIONS
// -------------------------------------------------------------------

void solve_stls(input in, bool verbose) {

  // Arrays for STLS solution
  float *xx = NULL; 
  float *phi = NULL;
  float *GG = NULL;
  float *GG_new = NULL;
  float *SS = NULL;
  float *SSHF = NULL;

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
  float iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    clock_t tic = clock();
    
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
  write_text(SS, GG, phi, SSHF, xx, in);
  write_guess(SS, GG, in); 
  if (verbose) printf("Done.\n");

  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);

 
}

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays(input in, float **xx, float **phi, 
		       float **GG, float **GG_new, 
		       float **SS, float **SSHF){

  *xx = (float*)malloc( sizeof(float) * in.nx);
  *phi = (float*)malloc( sizeof(float) * in.nx * in.nl);
  *SSHF = (float*)malloc( sizeof(float) * in.nx);
  *GG = (float*)malloc( sizeof(float) * in.nx);
  *GG_new = (float*)malloc( sizeof(float) * in.nx);
  *SS = (float*)malloc( sizeof(float) * in.nx);
  
}

void free_stls_arrays(float *xx, float *phi, float *GG, 
		      float *GG_new, float *SS,
		      float *SSHF){

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

void init_fixed_stls_arrays(input *in, float *xx, 
			    float *phi, float *SSHF, bool verbose){

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

  // Threads and blocks
  int threadLimitPerBlock = 1024;
  int numberOfThreads = in.nx * in.nl;
  int numberOfBlocks = (numberOfThreads/threadLimitPerBlock) + 1;
  dim3 grid(numberOfBlocks, 1, 1);
  dim3 block(threadLimitPerBlock, 1, 1);

  // Arrays for device
  float *d_phi = NULL;
  float *d_xx = NULL;
  cudaMalloc(&d_phi, sizeof(float) * in.nx * in.nl);
  cudaMalloc(&d_xx, sizeof(float) * in.nx);

  // Copy arrays to device
  cudaMemcpy(d_phi, psi, sizeof(float) * in.nx *in.nl, 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(d_xx, xx, sizeof(float) * in.nx, 
	     cudaMemcpyHostToDevice);  

  // Compute normalized ideal Lindhart density on device
  compute_phil <<< grid, block >>> (d_phi, d_xx, ll, in);
 
  // Copy solution from device 
  cudaMemcpy(phi, d_phi, sizeof(float)*in.nx*in.nl, cudaMemcpyDeviceToHost);

  // Free memory
  cudaFree(d_phi);
  cudaFree(d_xx);
 
}

__global__ void compute_phil(float *phil, float *xx,  int ll, input in) {

  int threadId = threadIdx.x + blockIdx.x * blockDim.x;

  if (threadId >= in.nx*in.nl) return;

  int ii = threadId % in.nx;
  int ll = threadId % in.nl;

  phil[threadId] = 0.0;
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

__device__ float phixl(float yy, float xx, int ll, input in) {

  float yy2 = yy*yy, xx2 = xx*xx, txy = 2*xx*yy, 
    tplT = 2*M_PI*ll*in.Theta, tplT2 = tplT*tplT;

  if (xx > 0.0) {
    return 1.0/(2*xx)*yy/(exp(yy2/in.Theta - in.mu) + 1.0)
      *log(((xx2+txy)*(xx2+txy) + tplT2)/((xx2-txy)*(xx2-txy) + tplT2));
  }
  else {
    return 0;
  }

}

__device__ float phix0(float yy, float xx, input in) {

  float yy2 = yy*yy, xx2 = xx*xx, xy = xx*yy;

  if (xx > 0.0){

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

  else{
    return (2.0/in.Theta)*yy2/(exp(yy2/in.Theta - in.mu)
                               + exp(-yy2/in.Theta + in.mu) + 2.0);
  }

}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf(float *SS, float *SSHF, float *GG, 
		 float *phi, float *xx, input in){

  float lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  float ff = 4*lambda*in.rs/M_PI;
  float xx2, BB, BB_tmp, BB_den, phixl;
  float tplT, Axl;

  for (int ii=0; ii<in.nx; ii++){

    if (xx[ii] > 0.0){
      xx2 = xx[ii]*xx[ii];
      BB = 0.0;
      
      for (int ll=0; ll<in.nl; ll++){
	tplT = 2*M_PI*ll*in.Theta;
	phixl = phi[idx2(ii,ll,in.nx)];
	Axl = (4.0/3.0)*xx2/(tplT*tplT + xx2*xx2);
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
 
float ssfHF(float yy, float xx, input in) {

  float yy2 = yy*yy, ypx = yy + xx, ymx = yy - xx;
 
  if (xx > 0.0){
    return -3.0*in.Theta/(4.0*xx)*yy/(exp(yy2/in.Theta - in.mu) + 1.0)
      *log((1 + exp(in.mu - ymx*ymx/in.Theta))
           /(1 + exp(in.mu - ypx*ypx/in.Theta)));
  }
  else {
    return -3.0/2.0*yy2/(1.0 + cosh(yy2/in.Theta - in.mu));
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

  if (xx > 0.0 && yy > 0.0){
    
    if (xx > yy){
      return -(3.0/4.0) * yy2 * (SS - 1.0)
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(xx - yy)));
    }
    else if (xx < yy) {
      return -(3.0/4.0) * yy2 * (SS - 1.0)
	* (1 + (xx2 - yy2)/(2*xx*yy)*log((xx + yy)/(yy - xx)));
    }
    else {
      return -(3.0/4.0) * yy2 * (SS - 1.0);
    }
    
  }
  else
    return 0;
  

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


// write text files for output
void write_text(float *SS, float *GG, float *phi, 
		float *SSHF, float *xx, input in){


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
    float lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
    float ff = 4*lambda*in.rs/M_PI;
    float sdr;
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

}


// write binary file to use as initial guess (or restart)
void write_guess(float *SS, float *GG, input in){

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
  fwrite(SS, sizeof(float), in.nx, fid);

  // Static local field correction
  fwrite(GG, sizeof(float), in.nx, fid);

  // Close binary file
  fclose(fid);

}


// read binary file to use as initial guess (or restart)
void read_guess(float *SS, float *GG, input in){

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
  fread(SS, sizeof(float), in_load.nx, fid);

  // Static structure factor in the Hartree-Fock approximation
  fread(GG, sizeof(float), in_load.nx, fid);

  // Close binary file
  fclose(fid);
	    
}
