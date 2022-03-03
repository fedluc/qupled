#ifndef STLS_H
#define STLS_H

#include "read_input.h"

// -------------------------------------------------------------------
// DATA STRUCTURE USED TO ALLOCATE THE STLS ARRAYS
// -------------------------------------------------------------------

typedef struct {

  double *xx;
  double *phi;
  double *SS;
  double *SSHF;
  double *GG;
  double *GG_new;
  
} stls_pointers;

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays_new(input in, stls_pointers *pp);

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
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE STLS SCHEME
// -------------------------------------------------------------------

void stls_iterations(double *SS, double *SSHF,
		     double *GG, double *GG_new,
		     double *phi, double *xx,
		     input in, bool verbose);

double stls_err(double *GG, double *GG_new, input in);

void stls_update(double *GG, double *GG_new, input in);

// -------------------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY AT FINITE TEMPERATURE
// -------------------------------------------------------------------------------------

void compute_idr(double *phi, double *xx, input in, bool verbose);

void compute_idr_one_frequency(double *phil, double *xx, int ll, input in);

double idr_partial_xl(double yy, void *pp);

double idr_partial_x0(double yy, void *pp);

// -------------------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY AT ZERO TEMPERATURE
// -------------------------------------------------------------------------------------

double idr_re_zero_temperature(double xx, double Omega);

double idr_im_zero_temperature(double xx, double Omega);

double idrp_re_zero_temperature(double xx, double Omega);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// ---------------------------------------------------------------------------

void compute_ssf_stls(double *SS, double *SSHF,
		      double *GG, double *phi, 
		      double *xx, input in);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC STRUCTURE FACTOR AT FINITE TEMPERATURE
// ---------------------------------------------------------------------------

void compute_ssf_stls_finite_temperature(double *SS, double *SSHF,
					 double *GG, double *phi, 
					 double *xx, input in);

void compute_ssf_HF_finite_temperature(double *SS,  double *xx, input in);

double ssf_HF_finite_temperature(double yy, void *pp);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC STRUCTURE FACTOR AT ZERO TEMPERATURE
// ---------------------------------------------------------------------------

void compute_ssf_stls_zero_temperature(double *SS, double *SSHF,
				       double *GG, double *xx,
				       input in);

void compute_ssf_HF_zero_temperature(double *SS,  double *xx, input in);

double ssf_stls_zero_temperature(double Omega, void* pp);

double ssf_plasmon(double xx, double GG, input in);

double drf_re_zero_temperature(double xx, void* pp);

double drfp_re_zero_temperature(double xx, double Omega,
				double GG, double rs);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc(double *GG, double *SS, 
		  double *xx, input in);

double slfc(double yy, void *pp);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_stls(double *SS, double *GG, double *phi, 
		       double *SSHF, double *xx, input in);

void write_text_ssf(double *SS, double *xx, input in);

void write_text_ssf_HF(double *SS, double *xx, input in);

void write_text_slfc(double *GG, double *xx, input in);

void write_text_sdr(double *GG, double *phi, double *xx, input in);

void write_text_idr(double *phi, input in);

void write_text_uint(double *SS, double *xx, input in);
  
void write_text_rdf(double *SS, double *xx, input in);

void write_guess_stls(double *SS, double *GG, input in);

void read_guess_stls(double *SS, double *GG, input in);

void check_guess_stls(int nx, double dx, double xmax,
		      input in, size_t it_read, size_t it_expected,
		      FILE *fid, bool check_grid, bool check_items,
		      bool check_eof);

#endif
