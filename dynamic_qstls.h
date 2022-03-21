#ifndef DYNAMIC_QSTLS_H
#define DYNAMIC_QSTLS_H

#include "read_input.h"

// -------------------------------------------------------------------
// CONSTANTS
// -------------------------------------------------------------------

// Number of data points for the integration of the second part
// auxiliary density response
#define ADR_NU 201

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

void compute_dynamic_adr_re_lev1(double *psi_re, double *WW,
				  double *SS, double *xx,
				  input in);

double adr_re_lev1_partial_xW(double ww, void* pp);

void compute_dynamic_adr_re_lev2(double *int_lev1, double WW,
				  double *ww, input in);

double adr_re_lev2_partial_xwW(double uu, void* pp);

void compute_dynamic_adr_re_lev3(double *int_lev2, double WW,
				  double ww, double *qq, double *uu,
				  input in);

double adr_re_lev3_partial_xwuW(double qq, void* pp);

double adr_re_lev3_partial_xwu0(double qq, void* pp);

// --------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE IMAGINARY PART OF THE AUXILIARY 
// DENSITY RESPONSE
// --------------------------------------------------------------------

void compute_dynamic_adr_im_lev1(double *psi_re, double *WW,
				  double *SS, double *xx,
				  input in);
double adr_im_lev1_partial_xW(double ww, void* pp);

void compute_dynamic_adr_im_lev2(double *psi_im_lev1, double WW,
				  double *ww, input in);

double adr_im_lev2_partial_xwW(double uu, void* pp);

double adr_im_lev2_partial_xw0(double uu, void* pp);

void compute_dynamic_adr_im_lev3(double *psi_im_lev2, double WW,
				  double ww, double *qq, double *uu,
				  input in);

double adr_im_lev3_partial_xwuW(double qq, void* pp);

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
