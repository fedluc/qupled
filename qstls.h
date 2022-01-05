#ifndef QSTLS_H
#define QSTLS_H

#include "solvers.h"

// ------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FIXED COMPONENT OF THE AUXILIARY RESPONSE
// ------------------------------------------------------------------------

void compute_adr_fixed(double *psi_fixed, double *xx, input in);

double adr_fixed_part1_partial_xwl(double qq, void* pp);

double adr_fixed_part1_partial_xw0(double qq, void* pp);

double adr_fixed_part2_partial_xwq0(double tt, void* pp);

double adr_fixed_part2_partial_xwql(double tt, void* pp);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE CHANGING COMPONENT OF THE AUXILIARY RESPONSE
// ---------------------------------------------------------------------------

void compute_adr(double *psi, double *psi_fixed, double *SS,
		 double *xx, input in);

double adr_partial(double ww, void* pp);


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_qstls(double *SS, double *SSHF, double *psi,
		       double *phi, double *xx, input in);


// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A THREE-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx3(int xx, int yy, int zz,
         int x_size, int y_size);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_qstls(double *SS, double *psi, double *phi, 
			double *SSHF, double *xx, input in);

void write_guess_qstls(double *SS, double *psi, input in);

void read_guess_qstls(double *SS, double *psi, input in);

void write_fixed_qstls(double *psi_fixed, input in);

void read_fixed_qstls(double *psi_fixed, input in);

#endif
