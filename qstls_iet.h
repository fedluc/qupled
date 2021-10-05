#ifndef QSTLS_IET_H
#define QSTLS_IET_H

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

void compute_psi_iet(double *psi_new, double *psi, double *psi_xlw_qstls,
                     double *phi, double *SS, double *bf, double *xx,
                     input in);

double psi_u_iet(double uu, void* pp);

double psi_w_iet(double ww, void* pp);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_qstls_iet(double *SS, double *SSHF, double *psi,
			   double *phi, double *bf, double *xx, input in);



// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_qstls_iet(double *SS, double *psi, double *phi,
                          double *SSHF, double *bf, double *xx,
                          input in);

#endif
