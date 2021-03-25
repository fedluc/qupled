#ifndef STLS_H
#define STLS_H

#include <stdbool.h>

// -------------------------------------------------------------------
// STRUCTURE TO STORE THE INPUT PARAMETERS
// -------------------------------------------------------------------

typedef struct {

  char *phi_file;
  char *ssf_file;
  float Theta;
  float rs;
  float dx;
  float err_min_iter;
  float a_mix;
  float mu_lo;
  float mu_hi;
  float mu;
  float xmax;
  int nl;
  int nx;
  int nIter;


} input;


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS EQUATIONS
// -------------------------------------------------------------------

void solve_stls(input in, bool verbose,
                float **xx_out, float **SS_out,
                float **SSHF_out, float **GG_out,
                float **GG_new_out, float **phi_out);

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays(input in, float **xx, float **phi,
		       float **GG, float **GG_new,
		       float **SS, float **SSHF);

void free_stls_arrays(float *xx, bool free_xx,
                      float *phi, bool free_phi,
                      float *GG, bool free_GG,
                      float *GG_new, bool free_GG_new,
                      float *SS, bool free_SS,
                      float *SSHF, bool free_SSHF);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

float compute_mu(input in);

double normalization_condition(double mu, void *pp);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------


void wave_vector_grid(float *xx, input in);

// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A TWO-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx2(int xx, int yy, int x_size);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY
// -------------------------------------------------------------------

void compute_phi(float *phi, float *xx, input in, bool verbose);

void compute_phil(float *phil, float *xx, int ll, input in);

float phixl(float yy, float xx, int ll, input in);

float phix0(float yy, float xx, input in);

void compute_ssfHF(float *SS,  float *xx, input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

float ssfHF(float yy, float xx, input in);


void compute_ssf(float *SS, float *SSHF,
                 float *GG, float *phi, 
		 float *xx, input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc(float *GG, float *SS, 
		  float *xx, input in);

float slfc(float yy, float xx, float SS);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

float compute_uex(float *SS, input in);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text(float *SS, float *GG, 
		float *xx, input in );

void write_bin(float *phi, float *SSHF, input in);

void read_bin(input *in, float **xx, float **phi,
              float **GG, float **GG_new,
              float **SS, float  **SSHF);

#endif
