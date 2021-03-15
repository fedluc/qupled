#ifndef STLS_HNC_H
#define STLS_HNC_H

#include "stls.h"

void solve_stls_hnc(input in, bool verbose);

void compute_slfc_hnc(float *GG_new, float *GG, 
		      float *SS, float *xx, input in);

void write_text_hnc(float *SS, float *GG,
		    float *xx, input in );


#endif
