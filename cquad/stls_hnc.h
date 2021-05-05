#ifndef STLS_HNC_H
#define STLS_HNC_H

#include "solvers.h"

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_hnc(double *GG_new, double *GG, double *SS,
                      double *bf, double *xx, input in);

double slfc_u(double uu, void* pp);

double slfc_w(double ww, void* pp);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE BRIDGE FUNCTION TERM
// -------------------------------------------------------------------

void compute_bf(double *bf, double *xx, input in);

void bf_hnc(double *bf, double *xx, input in);

void bf_ocp_ichimaru(double *bf, double *xx, input in);

void bf_ocp_2021(double *bf, double *xx, input in);

double rbfr(double rr, void *pp);

#endif
