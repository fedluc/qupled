#ifndef STLS_H
#define STLS_H

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

// ---------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE INITIAL GUESS
// ---------------------------------------------------------------------

void initial_guess_stls(double *xx, double *SS, double *SSHF,
			double *GG, double *GG_new, double *phi,
			input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS
// -------------------------------------------------------------------

double stls_err(double *GG, double *GG_new, input in);

void stls_update(double *GG, double *GG_new, input in);

// -------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY RESPONSE 
// -------------------------------------------------------------------------

void compute_idr(double *phi, double *xx, input in, bool verbose);

double idr_re_zero_temperature(double xx, double Omega);

double idr_im_zero_temperature(double xx, double Omega);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// ---------------------------------------------------------------------------

void compute_ssf_stls(double *SS, double *SSHF,
		      double *GG, double *phi, 
		      double *xx, input in);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE HARTREE-FOCK STATIC STRUCTURE FACTOR
// ---------------------------------------------------------------------------

void compute_ssf_HF(double *SS, double *xx, input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc(double *GG, double *SS, 
		  double *xx, input in);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_stls(double *SS, double *GG, double *phi, 
		       double *SSHF, double *xx, input in);

void write_text_ssf(double *SS, double *xx, input in);

void write_text_ssf_HF(double *SS, double *xx, input in);

void write_text_idr(double *phi, input in);

void write_text_uint(double *SS, double *xx, input in);
  
void write_text_rdf(double *SS, double *xx, input in);

void write_restart_stls(double *SS, double *GG, input in);

void read_restart_stls(double *SS, double *GG, input in);

#endif
