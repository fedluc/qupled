#ifndef DYNAMIC_QSTLS_H
#define DYNAMIC_QSTLS_H

#include "read_input.h"

// -------------------------------------------------------------------
// CONSTANTS
// -------------------------------------------------------------------

// Number of data points for the integration of the second part
// of the imaginary auxiliary density response
#define ADR_IM_NU 201

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC PROPERTIES OF THE QSTLS
// SCHEME
// -------------------------------------------------------------------

void compute_dynamic_qstls(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_qstls_arrays(input in, double **psi_re, 
				double **psi_im);
void free_dynamic_qstls_arrays(double *psi_re, double *psi_im, double *SS,
			       double *xx);

// -------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE STATIC STRUCTURE FACTOR (FROM FILE)
// -------------------------------------------------------------------

void get_ssf(double **SS, double **xx, input *in);

// --------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE REAL PART OF THE AUXILIARY DENSITY
// RESPONSE
// --------------------------------------------------------------------

void compute_dynamic_adr(double *psi_re, double *psi_im,
			 double *WW, double *SS,
			 double *xx, input in);

void compute_dynamic_adr_re_part1(double *psi_re, double *WW,
				  double *SS, double *xx,
				  input in);

double adr_re_part1_partial_xW(double ww, void* pp);

void compute_dynamic_adr_re_part2(double *psi_re_part1, double WW,
				  double *ww, input in);

double adr_re_part2_partial_xwW(double qq, void* pp);

double adr_re_part2_partial_xw0(double qq, void* pp);

void compute_dynamic_adr_re_part3(double *psi_re_part2, double WW,
				  double ww, double *qq, input in);

double adr_re_part3_partial_xwqW(double tt, void* pp);

double adr_re_part3_partial_xwq0(double tt, void* pp);

// --------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE IMAGINARY PART OF THE AUXILIARY 
// DENSITY RESPONSE
// --------------------------------------------------------------------

void compute_dynamic_adr_im_part1(double *psi_re, double *WW,
				  double *SS, double *xx,
				  input in);
double adr_im_part1_partial_xW(double ww, void* pp);

void compute_dynamic_adr_im_part2(double *psi_im_part1, double WW,
				  double *ww, input in);

double adr_im_part2_partial_xwW(double uu, void* pp);

void compute_dynamic_adr_im_part3(double *psi_im_part2, double WW,
				  double ww, double *qq, double *uu,
				  input in);

double adr_im_part3_partial_xwuW(double qq, void* pp);

// ---------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC STRUCTURE FACTOR
// ---------------------------------------------------------------------

void compute_dsf_qstls(double *SSn, double *phi_re, double *phi_im,
		       double *psi_re, double *psi_im,
		       double *WW, input in);


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_dynamic_qstls(double *SSn, double *WW, double *psi_re,
			      double *psi_im, input in);


#endif
