#ifndef QSTLS_H
#define QSTLS_H

#include "stls.h"

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE QSTLS EQUATIONS
// -------------------------------------------------------------------

void solve_qstls(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_qstls_arrays(input in, double **psi, double **SS_new);

void free_qstls_arrays(double *xx, double *phi,
		       double *psi, double *SS,
		       double *SS_new, double *SSHF);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE AUXILLIARY RESPONSE
// -------------------------------------------------------------------

void compute_psi(double *psi, double *xx, 
		 double *SS,  input in);

void compute_psil(double *psil, double *xx,  double *SS, 
		  int ll, input in);

double psi_u(double uu, double qq, double ww,
             double xx, int ll, input in);

double psi_q(double qq, int ll, input in);

double psi_w(double ww, double SS);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_qstls_ssf(double *SS, double *SSHF, double *phi,
                       double *psi, double *xx, input in);

// -------------------------------------------------------------------
// FUNCTION FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_qstls(double *SS, double *phi, double *psi, double *xx, input in );

#endif
