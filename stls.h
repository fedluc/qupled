#ifndef STLS_H
#define STLS_H

#include "solvers.h"
#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays(input in, double **xx, double **phi,
		       double **GG, double **GG_new,
		       double **SS, double **SSHF);

void free_stls_arrays(double *xx, double *phi, double *GG,
                      double *GG_new, double *SS,
                      double *SSHF);


// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_stls_arrays(input *in, double *xx, 
			    double *phi, double *SSHF, bool verbose);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void wave_vector_grid(double *xx, input *in);

// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A TWO-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx2(int xx, int yy, int x_size);


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY
// -------------------------------------------------------------------

void compute_idr(double *phi, double *xx, input in, bool verbose);

void compute_idr_one_frequency(double *phil, double *xx, int ll, input in);

double idr_partial_xl(double yy, void *pp);

double idr_partial_x0(double yy, void *pp);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_stls(double *SS, double *SSHF,
		      double *GG, double *phi, 
		      double *xx, input in);

void compute_ssf_HF(double *SS,  double *xx, input in);

double ssf_HF(double yy, void *pp);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc(double *GG, double *SS, 
		  double *xx, input in);

double slfc(double yy, void *pp);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

double compute_internal_energy(double *SS, double *xx, input in);

double uex(double yy, void *pp);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_stls(double *SS, double *GG, double *phi, 
		       double *SSHF, double *xx, input in);

void write_guess_stls(double *SS, double *GG, input in);

void read_guess_stls(double *SS, double *GG, input in);

void check_guess_stls(int nx, double dx, double xmax, input in,
		      size_t it_read, size_t it_expected, FILE *fid,
		      bool check_grid, bool check_items, bool check_eof);

#endif
