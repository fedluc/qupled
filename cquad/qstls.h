#ifndef QSTLS_H
#define QSTLS_H

#include "solvers.h"

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE ...
// -------------------------------------------------------------------

void compute_psi_xlw(double *psi_xlw, double *xx, input in);

double psi_x0w_t(double uu, void* pp);
double psi_x0w_q(double uu, void* pp);
double psi_xlw_t(double uu, void* pp);
double psi_xlw_q(double uu, void* pp);

#endif
