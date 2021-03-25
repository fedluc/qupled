#ifndef STLS_HNC_H
#define STLS_HNC_H

#include "stls.h"

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS-HNC EQUATIONS
// -------------------------------------------------------------------

void solve_stls_hnc(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_hnc(float *GG_new, float *GG, 
		      float *SS, float *xx, input in);

// -------------------------------------------------------------------
// FUNCTION FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_hnc(float *SS, float *GG,
		    float *xx, input in );


#endif
