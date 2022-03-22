#ifndef DYNAMIC_QSTLS_IET_H
#define DYNAMIC_QSTLS_IET_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC PROPERTIES OF THE QSTLS
// SCHEME
// -------------------------------------------------------------------

void compute_dynamic_qstls_iet(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_qstls_iet_arrays(input in, double **phi_re_2D,
				    double **phi_im_2D,
				    double **psi_re_2D,
				    double **psi_im_2D,
				    double **phi_re_1D,
				    double **phi_im_1D,
				    double **psi_re_1D,
				    double **psi_im_1D);

void free_dynamic_qstls_iet_arrays(double *phi_re_2D,
				   double *phi_im_2D,
				   double *psi_re_2D,
				   double *psi_im_2D,
				   double *phi_re_1D,
				   double *phi_im_1D,
				   double *psi_re_1D,
				   double *psi_im_1D);

// -------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE STATIC STRUCTURE FACTOR (FROM FILE)
// -------------------------------------------------------------------

void get_ssf(double **SS, double **xx, input *in);

// --------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE IDEAL DENSITY RESPONSE
// --------------------------------------------------------------------

void compute_dynamic_idr_iet(double *phi_re, double *phi_im,
			     double *WW, double *xx, input in);

// -------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE BRIDGE FUNCTION
// -------------------------------------------------------------------

void get_bf(double **bf, double *xx, input in);

// --------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE AUXILIARY DENSITY RESPONSE
// --------------------------------------------------------------------

void compute_dynamic_adr_iet(double *psi_re, double *psi_im,
			     double *phi_re, double *phi_im,
			     double *WW, double *SS,
			     double *bf, double *xx,
			     input in);

// --------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE REAL PART OF THE AUXILIARY DENSITY
// RESPONSE
// --------------------------------------------------------------------

void compute_dynamic_adr_iet_re(double *psi_re, double *phi_re,
				double *WW, double *SS,
				double *bf, double *xx,
				input in);

double adr_iet_err(double *psi_re, double *psi_re_new, input in);

void adr_iet_update(double *psi_re, double *psi_re_new, input in);

void compute_dynamic_adr_iet_re_lev1(double *psi_re_new, double *psi_re,
				     double *psi_re_fixed, double *phi_re,
				     double *WW, double *SS, double *bf,
				     double *xx, input in);

void compute_dynamic_adr_iet_re_lev1_1(double *int_lev1_1, double *psi_re,
				       double *phi_re, double *SS,
				       double *bf, input in);

void compute_dynamic_adr_iet_re_lev1_2(double *int_lev1_2,
				       double *psi_re_fixed,
				       int ii, int jj,
				       input in, bool read);
  
double adr_iet_re_lev1_partial_xW(double ww, void* pp);

void compute_dynamic_adr_iet_re_lev2(double *int_lev1, double WW,
				     double xx, double *SS,
				     double *ww, input in);

double adr_iet_re_lev2_partial_xwW(double uu, void* pp);

void compute_dynamic_adr_iet_re_lev3(double *int_lev2, double WW,
				     double xx, double ww,
				     double *uu, input in);

double adr_iet_re_lev3_partial_xwuW(double qq, void* pp);

double adr_iet_re_lev3_partial_xwu0(double qq, void* pp);


// --------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE IMAGINARY PART OF THE AUXILIARY
// DENSITY RESPONSE
// --------------------------------------------------------------------

void compute_dynamic_adr_iet_im_lev1(double *psi_im, double *psi_re,
				     double *phi_re, double *WW,
				     double *SS, double *bf,
				     double *xx, input in);

void compute_dynamic_adr_iet_im_lev1_1(double *int_lev1_1, double *psi_re,
				       double *phi_re, double *SS,
				       double *bf, input in);

double adr_iet_im_lev1_partial_xW(double ww, void* pp);

void compute_dynamic_adr_iet_im_lev2(double *int_lev1, double WW,
				     double xx, double *SS,
				     double *ww, input in);

double adr_iet_im_lev2_partial_xwW(double uu, void* pp);

double adr_iet_im_lev2_partial_xw0(double uu, void* pp);

void compute_dynamic_adr_iet_im_lev3(double *int_lev2, double WW,
				     double xx, double ww,
				     double *uu, input in);

double adr_iet_im_lev3_partial_xwuW(double qq, void* pp);

// ---------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC STRUCTURE FACTOR
// ---------------------------------------------------------------------

void compute_dsf_qstls_iet(double *SSn, double *phi_re, double *phi_im,
			   double *psi_re, double *psi_im,
			   double *WW ,input in);


/* // ------------------------------------------------------------------- */
/* // FUNCTIONS FOR OUTPUT AND INPUT */
/* // ------------------------------------------------------------------- */

/* void write_text_dynamic_qstls(double *SSn, double *WW, double *psi_re, */
/* 			      double *psi_im, input in); */


#endif
