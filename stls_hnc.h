#ifndef STLS_HNC_H
#define STLS_HNC_H

#include "stls.h"

void solve_stls_hnc(input in, bool verbose);

void compute_slfc_hnc(double *GG, double *SS, 
		  double *xx, input in);

void write_text_hnc(double *SS, double *GG,
		    double *xx, input in );


double slfc_u(double yy, void* pp);
double slfc_w(double yy, void* pp);

#endif
