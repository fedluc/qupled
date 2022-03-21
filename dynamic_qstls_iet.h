#ifndef DYNAMIC_QSTLS_IET_H
#define DYNAMIC_QSTLS_IET_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC PROPERTIES OF THE QSTLS
// SCHEME
// -------------------------------------------------------------------

void compute_dynamic_qstls_iet(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE STATIC STRUCTURE FACTOR (FROM FILE)
// -------------------------------------------------------------------

void get_ssf(double **SS, double **xx, input *in);

// -------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE BRIDGE FUNCTION
// -------------------------------------------------------------------

void get_bf(double **bf, double *xx, input in);

// --------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE REAL PART OF THE AUXILIARY DENSITY
// RESPONSE
// --------------------------------------------------------------------

void compute_dynamic_adr_iet(double *psi_re, double *psi_im,
			     double *WW, double *SS,
			     double *bf, double *xx,
			     input in);

void compute_dynamic_adr_re_iet(double *psi_re, double *WW,
				double *SS, double *bf,
				double *xx, input in);

void compute_dynamic_adr_iet_re_lev1(double *psi_re, double psi_re_old,
				      double *WW, double *SS,
				      double *bf, double *xx,
				      input in);

double adr_iet_re_lev1_partial_xW(double ww, void* pp);

void compute_dynamic_adr_iet_re_lev2(double *int_lev1, double WW,
				      double *ww, input in);

double adr_iet_re_lev2_partial_xwW(double uu, void* pp);

void compute_dynamic_adr_iet_re_lev3(double *int_lev2, double WW,
				      double ww, double *qq, double *uu,
				      input in);

double adr_iet_re_lev3_partial_xwuW(double qq, void* pp);

double adr_iet_re_lev3_partial_xwu0(double qq, void* pp);

/* // -------------------------------------------------------------------- */
/* // FUNCTIONS USED TO COMPUTE THE IMAGINARY PART OF THE AUXILIARY  */
/* // DENSITY RESPONSE */
/* // -------------------------------------------------------------------- */

/* void compute_dynamic_adr_im_lev1(double *psi_re, double *WW, */
/* 				  double *SS, double *xx, */
/* 				  input in); */
/* double adr_im_lev1_partial_xW(double ww, void* pp); */

/* void compute_dynamic_adr_im_lev2(double *psi_im_lev1, double WW, */
/* 				  double *ww, input in); */

/* double adr_im_lev2_partial_xwW(double uu, void* pp); */

/* double adr_im_lev2_partial_xw0(double uu, void* pp); */

/* void compute_dynamic_adr_im_lev3(double *psi_im_lev2, double WW, */
/* 				  double ww, double *qq, double *uu, */
/* 				  input in); */

/* double adr_im_lev3_partial_xwuW(double qq, void* pp); */

/* // --------------------------------------------------------------------- */
/* // FUNCTION USED TO COMPUTE THE DYNAMIC STRUCTURE FACTOR */
/* // --------------------------------------------------------------------- */

/* void compute_dsf_qstls(double *SSn, double *phi_re, double *phi_im, */
/* 		       double *psi_re, double *psi_im, */
/* 		       double *WW, input in); */


/* // ------------------------------------------------------------------- */
/* // FUNCTIONS FOR OUTPUT AND INPUT */
/* // ------------------------------------------------------------------- */

/* void write_text_dynamic_qstls(double *SSn, double *WW, double *psi_re, */
/* 			      double *psi_im, input in); */


#endif
