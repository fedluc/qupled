#ifndef STLS_IET_H
#define STLS_IET_H

#include "solvers.h"

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

void write_bridge_function(double *bf, double *xx, input in);

#endif
