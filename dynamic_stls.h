#ifndef DYNAMIC_STLS_H
#define DYNAMIC_STLS_H

#include "read_input.h"


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC PROPERTIES OF THE CLASSICAL
// SCHEMES (STLS, VS-STLS AND STLS-IET)
// -------------------------------------------------------------------

void compute_dynamic_stls(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE SIZE OF THE FREQUENCY GRID
// -------------------------------------------------------------------

void get_frequency_grid_size(input *in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_stls_arrays(input in, double **WW, double **phi_re, 
			       double **phi_im, double **SSn);

void free_dynamic_stls_arrays(double *WW, double *phi_re, 
			      double *phi_im, double *SSn);

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_dynamic_stls_arrays(input *in, double *WW,
				    bool verbose);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE FREQUENCY GRID
// ------------------------------------------------------------------

void frequency_grid(double *WW, input *in);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE IDEAL DENSITY RESPONSE
// ------------------------------------------------------------------

void compute_dynamic_idr(double *phi_re, double *phi_im,
			 double *WW, input in);

void compute_dynamic_idr_re(double *phi_re, double *WW,
			    input in);

void compute_dynamic_idr_im(double *phi_im, double *WW,
			    input in);

double idr_re_partial_xw(double yy, void *pp);

double idr_re_partial_x0(double yy, void *pp);

double idr_im_partial_xw(double yy, void *pp);

double idr_im_partial_x0(double yy, void *pp);

// ---------------------------------------------------------------------
// FUNCTION USED TO OBTAIN THE STATIC LOCAL FIELD CORRECTION (FROM FILE)
// ---------------------------------------------------------------------

void get_slfc(double *GG, input in);

// ---------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE DYNAMIC STRUCTURE FACTOR
// ---------------------------------------------------------------------

void compute_dsf(double *SSn, double *phi_re, double *phi_im,
		 double GG, double *WW, input in);

// ---------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// ---------------------------------------------------------------------

void write_text_dynamic_stls(double *SSn, double *WW, input in);

void write_text_dsf(double *SSn, double *WW, input in);


// --------------------------------------------------------------------
// QSTLS
// --------------------------------------------------------------------

void compute_dynamic_qstls(input in, bool verbose);
void alloc_dynamic_qstls_arrays(input in, double **psi_re, 
				double **psi_im);
void free_dynamic_qstls_arrays(double *psi_re, double *psi_im, double *SS,
			       double *xx);
void get_ssf(double **SS, double **xx, input *in);
void compute_dynamic_adr(double *psi_re, double *psi_im,
			 double *WW, double *SS,
			 double *xx, input in);
void compute_dynamic_adr_re(double *psi_re, double *WW,
			    double *SS, double *xx,
			    input in);
double adr_re_part1_partial_xW(double ww, void* pp);
double adr_re_part2_partial_xwW(double qq, void* pp);
double adr_re_part2_partial_xw0(double qq, void* pp);
double adr_re_part3_partial_xwqW(double tt, void* pp);
double adr_re_part3_partial_xwq0(double tt, void* pp);

#endif
