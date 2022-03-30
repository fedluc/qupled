#ifndef DYNAMIC_QSTLS_H
#define DYNAMIC_QSTLS_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_qstls_arrays(input in, double **psi_re,
				double **psi_im);

void alloc_dynamic_qstls_2Darrays(input in, double **psi_re,
				  double **psi_im);

void free_dynamic_qstls_arrays(double *psi_re, double *psi_im,
			       double *SS, double *xx);

void free_dynamic_qstls_2Darrays(double *psi_re, double *psi_im);

// -------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE STATIC STRUCTURE FACTOR (FROM FILE)
// -------------------------------------------------------------------

void get_ssf(double **SS, double **xx, input *in);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_dynamic_qstls(double *psi_re, double *psi_im,
			      double *WW, input in);

void write_bin_dynamic_adr_2D(double *psi_re, double *psi_im,
			      input in);

void read_bin_dynamic_adr_2D(double *psi_re, double *psi_im,
			     input in);

#endif
