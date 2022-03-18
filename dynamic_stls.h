#ifndef DYNAMIC_STLS_H
#define DYNAMIC_STLS_H

#include "read_input.h"

// -------------------------------------------------------------------
// CONSTANTS
// -------------------------------------------------------------------

// Number of data points for imaginary time (varies between 0 and 1)
#define ISF_NTAU 100

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
// FUNCTION USED TO COMPUTE THE INTERMEDIATE SCATTERING FUNCTION
// ---------------------------------------------------------------------

void compute_isf(double *FF, double *tt, double *SSn,
		 double *WW, input in);
 
double isf(double WW, void *pp);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// ---------------------------------------------------------------------

void write_text_dynamic_stls(double *SSn, double *WW, input in);

void write_text_dsf(double *SSn, double *WW, input in);


void write_text_isf(double *SSn, double *ww, input in);

#endif
