#ifndef STLS_HNC_H
#define STLS_HNC_H

#include "solvers.h"

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_hnc(float *GG_new, float *GG, float *SS,
                      float *bf, float *xx, input in);

void compute_bf(float *bf, float *xx, input in, bool iet);


#endif
