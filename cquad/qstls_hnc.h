#ifndef QSTLS_HNC_H
#define QSTLS_HNC_H

#include "solvers.h"

// ------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FIXED COMPONENT OF THE AUXILIARY RESPONSE
// ------------------------------------------------------------------------

void compute_psi_xluw(double *xx, input in);

double psi_x0uw_y(double yy, void* pp);

double psi_xluw_y(double yy, void* pp);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE CHANGING COMPONENT OF THE AUXILIARY RESPONSE
// ---------------------------------------------------------------------------

void compute_psi_hnc(double *psi_new, double *psi, double *psi_xlw_qstls,
                     double *phi, double *SS, double *bf, double *xx,
                     input in);

double psi_u_hnc(double uu, void* pp);

double psi_w_hnc(double ww, void* pp);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

#endif
