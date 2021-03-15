extern "C" {
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include "stls.h"
#include "qstls.h"
#include "qstls_gpu.h"
}

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS-HNC EQUATIONS
// -------------------------------------------------------------------

void solve_qstls(input in, bool verbose) {

  // Arrays for STLS solution
  float *xx = NULL; 
  float *phi = NULL;
  float *psi = NULL;
  float *SS = NULL;
  float *SS_new = NULL;
  float *SSHF = NULL;
  float *GG = NULL;
  float *GG_new = NULL;

  // Solve STLS equation for initial guess
  if (verbose) printf("Solution of classical STLS for initial guess:\n");
  float a_mix_hold = in.a_mix;
  in.a_mix = 0.1;
  solve_stls(in, false, &xx, &SS, &SSHF, &GG, &GG_new, &phi);
  free(GG);
  free(GG_new);
  in.a_mix = a_mix_hold;
  if (verbose) printf("Done.\n");

  // Allocate arrays for qstls calculation
  alloc_qstls_arrays(in, &psi, &SS_new);  
  
  // SSF and SLFC via iterative procedure
  if (verbose) printf("SSF and SLFC calculation...\n");
  float iter_err = 1.0;
  int iter_counter = 0;
  while (iter_counter < in.nIter && iter_err > in.err_min_iter ) {
    
    // Start timing
    clock_t tic = clock();
    
    // Update auxilliary response 
    compute_psi(psi, xx, SS, in);
    
    // Update SSF
    compute_qstls_ssf(SS_new, SSHF, phi, psi, xx, in);

    // Update diagnostic
    iter_err = 0.0;
    iter_counter++;
    for (int ii=0; ii<in.nx; ii++) {
      iter_err += (SS_new[ii] - SS[ii]) * (SS_new[ii] - SS[ii]);
      SS[ii] = in.a_mix*SS_new[ii] + (1-in.a_mix)*SS[ii];
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
  write_text_qstls(SS, psi, xx, in);
  if (verbose) printf("Done.\n");

  // Output to variable or free memory
  free_qstls_arrays(xx, phi, psi, SS, SS_new, SSHF);
 
 
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_qstls_arrays(input in, float **psi, float **SS_new){

  *psi = (float*)malloc( sizeof(float) * in.nx * in.nl);
  *SS_new = (float*)malloc( sizeof(float) * in.nx);

}

void free_qstls_arrays(float *xx, float *phi,
                      float *psi, float *SS,
                      float *SS_new, float *SSHF){

  free(xx);
  free(phi);
  free(psi);
  free(SSHF);
  free(SS);
  free(SS_new);

}


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE AUXILLIARY RESPONSE
// -------------------------------------------------------------------

void compute_psi(float *psi, float *xx,  float *SS, input in) {

  // Threads and blocks
  int threadLimitPerBlock = 1024;
  int numberOfThreads = in.nx * in.nl;
  int numberOfBlocks = (numberOfThreads/threadLimitPerBlock) + 1;
  dim3 grid(numberOfBlocks, 1, 1);
  dim3 block(threadLimitPerBlock, 1, 1);

  // Arrays for device
  float *d_psi = NULL;
  float *d_xx = NULL;
  float *d_SS = NULL;
  cudaMalloc(&d_psi, sizeof(float) * in.nx * in.nl);
  cudaMalloc(&d_xx, sizeof(float) * in.nx);
  cudaMalloc(&d_SS, sizeof(float) * in.nx);

  // Copy arrays to device
  cudaMemcpy(d_psi, psi, sizeof(float) * in.nx *in.nl, 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(d_xx, xx, sizeof(float) * in.nx, 
	     cudaMemcpyHostToDevice);  
  cudaMemcpy(d_SS, SS, sizeof(float) * in.nx, 
	     cudaMemcpyHostToDevice);

  // Launch kernel on device
  compute_psil <<< grid, block >>> (d_psi, d_xx, d_SS, in);
 
  // Copy solution from device 
  cudaMemcpy(psi, d_psi, sizeof(float)*in.nx*in.nl, cudaMemcpyDeviceToHost);

  // Free memory
  cudaFree(d_psi);
  cudaFree(d_xx);
  cudaFree(d_SS);

}

__global__ void compute_psil(float *psi, float *xx, float *SS, input in) {

  int threadId = threadIdx.x + blockIdx.x * blockDim.x;

  if (threadId >= in.nx*in.nl) return;

  float nu = (int)floor(2.0/in.dx);
  float fw, fq;
  int ii = threadId % in.nx;
  int ll = threadId / in.nl;

  psi[threadId] = 0.0;
  for (int jj=0; jj<in.nx; jj++){
      
    fw = 0.0;

    for (int kk=0; kk<in.nx; kk++){

      fq = 0.0;
	
      for (int mm=0; mm<nu; mm++){
  	fq += psi_u(-1.0 + mm*in.dx, xx[kk], xx[jj],
  		    xx[ii], ll, in);
      }

      fw += fq*in.dx*psi_q(xx[kk], ll, in);
		  
    }

    psi[threadId] += fw*in.dx*psi_w(xx[ii], SS[ii]);
    
  }
  
  if (ll == 0) psi[threadId] *= -(3.0/8.0)*in.dx;
  else psi[threadId] *= -3.0/(4.0*in.Theta)*in.dx;
  
}



__device__ float psi_u(float uu, float qq, float ww,
	     float xx, int ll, input in){

  float xx2 = xx*xx, qq2 = qq*qq, ww2 = ww*ww, 
    txq = 2.0*xx*qq, tplT = 2*M_PI*ll*in.Theta, 
    tplT2 = tplT*tplT, xwu = xx*ww*uu, 
    tt = xx2 - xwu, tt2 = tt*tt, 
    fact = ww*xx/(ww2 + tt - xwu);
  float logarg;
    
  if (ll == 0){

    if (tt == txq || tt == -txq) {
      return fact * 2.0 * qq*tt/xx;
    }
    else {
      
      logarg = (tt + txq)/(tt - txq);
      if (logarg < 0) logarg = -logarg;
      return fact * ((qq2 - tt2/(4.0*xx2))* 
      		     log(logarg) + qq*tt/xx);
    }

  }

  else {

    return fact * log(((txq + tt)*(txq + tt) + tplT2)
		      /((txq - tt)*(txq - tt) + tplT2));

  }

}

__device__ float psi_q(float qq, int ll, input in){

  float qq2 = qq*qq;

  if (ll == 0){

    return qq/(exp(qq2/in.Theta - in.mu) 
		+ exp(-qq2/in.Theta + in.mu) + 2.0);

  }
  else {

    return qq/(exp(qq2/in.Theta - in.mu) +  1.0);

  }

}

__device__ float psi_w(float ww, float SS){

  return ww * (SS - 1);
  
}

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_qstls_ssf(float *SS, float *SSHF, float *phi,
		       float *psi, float *xx, input in){

  float lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);
  float pilambda = M_PI*lambda;
  float ff = 4*lambda*lambda*in.rs;
  float ff3_2T = 3.0*in.Theta*ff/2.0;
  float xx2, BB, BB_tmp, phixl, psixl;
  for (int ii=0; ii<in.nx; ii++){

    xx2 = xx[ii]*xx[ii];
    BB = 0.0;

    for (int ll=0; ll<in.nl; ll++){
      phixl = phi[idx2(ii,ll,in.nx)];
      psixl = psi[idx2(ii,ll,in.nx)];
      BB_tmp = phixl*(phixl - psixl)/
	(pilambda*xx2 + ff*(phixl-psixl));
      if (ll>0) BB_tmp *= 2.0;
      BB += BB_tmp;
    }
    
    SS[ii] = SSHF[ii] - ff3_2T*BB;

  }

}

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------


// write text files with SSF and SLFC
void write_text_qstls(float *SS, float *psi, float *xx, input in){


    FILE* fid;
    
    // Output for SSF
    fid = fopen("ssf_QSTLS.dat", "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the static structure factor");
        exit(EXIT_FAILURE);
    }
    for (int ii=0; ii<in.nx; ii++)
    {
        fprintf(fid, "%.8e %.8e\n", xx[ii], SS[ii]);
    }
    fclose(fid);

    // Output for auxilliary response
    fid = fopen("psi_QSTLS.dat", "w");
    if (fid == NULL) {
        perror("Error while creating the output file for the auxilliary response");
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

