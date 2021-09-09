#ifndef STLS_IET_H
#define STLS_IET_H

#include "solvers.h"

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_iet(double *GG_new, double *GG, double *SS,
                      double *bf, double *xx, input in);

double slfc_u(double uu, void* pp);

double slfc_w(double ww, void* pp);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE BRIDGE FUNCTION TERM
// -------------------------------------------------------------------

void compute_bf(double *bf, double *xx, input in);

void bf_hnc(double *bf, double *xx, input in);

void bf_ocp_ioi(double *bf, double *xx, input in);

void bf_ocp_lct(double *bf, double *xx, input in);

void bf_rescaled_ocp_lct(double *bf, double *xx, input in);

double rbfr(double rr, void *pp);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_bf_static(double *bf, double *xx, input in);

#endif
