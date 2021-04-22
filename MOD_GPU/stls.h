#ifndef STLS_H
#define STLS_H

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays(input in, float **xx, float **phi,
		       float **GG, float **GG_new,
		       float **SS, float **SSHF);

void free_stls_arrays(float *xx, float *phi, float *GG,
                      float *GG_new, float *SS,
                      float *SSHF);


// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_stls_arrays(input *in, float *xx, 
			    float *phi, float *SSHF, bool verbose);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

float compute_mu(input in);

float normalization_condition(float mu, void *pp);

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

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssfHF(float *SS,  float *xx, input in);

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

void write_text(float *SS, float *GG, float *phi, 
		float *SSHF, float *xx, input in);

void write_guess(float *SS, float *GG, input in);

void read_guess(float *SS, float *GG, input in);

#endif
