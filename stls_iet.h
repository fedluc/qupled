#ifndef STLS_IET_H
#define STLS_IET_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_iet_arrays(input in, double **bf);

void free_stls_iet_arrays(double *bf);

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_stls_iet_arrays(double *bf, double *xx, input in);

// -------------------------------------------------------------------
// FUNCTION USED TO PERFORM THE ITERATIONS FOR THE STLS SCHEME
// -------------------------------------------------------------------

void stls_iet_iterations(double *SS, double *SSHF,
			 double *GG, double *GG_new,
			 double *phi, double *bf,
			 double *xx, input in,
			 bool verbose);
  
// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_iet(double *GG_new, double *GG, double *SS,
                      double *bf, double *xx, input in);

double slfc_partial_part1(double uu, void* pp);

double slfc_partial_part2(double ww, void* pp);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE BRIDGE FUNCTION TERM
// -------------------------------------------------------------------

void compute_bridge_function(double *bf, double *xx, input in);

void bridge_function_hnc(double *bf, double *xx, input in);

void bridge_function_ocp_ioi(double *bf, double *xx, input in);

void bridge_function_ocp_lct(double *bf, double *xx, input in);

void bridge_function_rescaled_ocp_lct(double *bf, double *xx, input in);

double rbfr(double rr, void *pp);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_stls_iet(double *SS, double *GG, double *phi, 
			 double *SSHF, double *xx, double *bf,
			 input in);
  
void write_text_bf(double *bf, double *xx, input in);

#endif
