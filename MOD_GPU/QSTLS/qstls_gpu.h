#ifndef QSTLS_GPU_H
#define QSTLS_GPU_H

#include "stls.h"

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_qstls_arrays(input in, float **psi, float **SS_new);

void free_qstls_arrays(float *xx, float *phi,
		       float *psi, float *SS,
		       float *SS_new, float *SSHF);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE AUXILLIARY RESPONSE
// -------------------------------------------------------------------

void compute_psi(float *psi, float *xx, 
		 float *SS,  input in);

__global__ void compute_psil(float *psil, float *xx,  float *SS, 
			     input in);

__device__ float psi_u(float uu, float qq, float ww,
			float xx, int ll, input in);

__device__ float psi_q(float qq, int ll, input in);

__device__ float psi_w(float ww, float SS);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_qstls_ssf(float *SS, float *SSHF, float *phi,
                       float *psi, float *xx, input in);

// -------------------------------------------------------------------
// FUNCTION FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_qstls(float *SS, float *phi, float *psi, float *xx, input in );

#endif
