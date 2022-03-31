#ifndef DYNAMIC_STLS_H
#define DYNAMIC_STLS_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE SIZE OF THE FREQUENCY GRID
// -------------------------------------------------------------------

void get_frequency_grid_size(input *in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_dynamic_stls_arrays(input in, double **WW, double **phi_re,
			       double **phi_im, double **SSn);

void alloc_dynamic_stls_2Darrays(input in, double **phi_re,
				 double **phi_im);

void free_dynamic_stls_arrays(double *WW, double *phi_re,
			      double *phi_im, double *SSn,
			      double *xx);

void free_dynamic_stls_2Darrays(double *phi_re,
				double *phi_im);

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_dynamic_stls_arrays(input *in, double *WW,
				    bool verbose);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE IDEAL DENSITY RESPONSE
// ------------------------------------------------------------------

void compute_dynamic_idr(double *phi_re, double *phi_im,
			 double *WW, double *xx, input in);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// ---------------------------------------------------------------------

void write_text_dynamic_stls(double *SSn, double *phi_re,
			     double *phi_im, double *WW,
			     input in);

void check_bin_dynamic(double dx, int nx, double xmax,
		       double dW, int dyn_nW, double Wmax,
		       double Wmin, double Theta, double rs,
		       input in, size_t it_read,
		       size_t it_expected, FILE *fid,
		       bool check_grid, bool check_items,
		       bool check_eof);
#endif
