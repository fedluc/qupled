#ifndef STLS_H
#define STLS_H

#include <stdbool.h>

// -------------------------------------------------------------------
// STRUCTURE TO STORE THE INPUT PARAMETERS
// -------------------------------------------------------------------

typedef struct {

  char *guess_file;
  char *theory;
  double Theta;
  double rs;
  double dx;
  double err_min_iter;
  double a_mix;
  double mu_lo;
  double mu_hi;
  double mu;
  double xmax;
  int nl;
  int nx;
  int nIter;


} input;

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS EQUATIONS
// -------------------------------------------------------------------

void solve_stls(input in, bool verbose,
                double **xx_out, double **SS_out,
                double **SSHF_out, double **GG_out,
                double **GG_new_out, double **phi_out);

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_stls_arrays(input in, double **xx, double **phi,
		       double **GG, double **GG_new,
		       double **SS, double **SSHF);

void free_stls_arrays(double *xx, bool free_xx,
                      double *phi, bool free_phi,
                      double *GG, bool free_GG,
                      double *GG_new, bool free_GG_new,
                      double *SS, bool free_SS,
                      double *SSHF, bool free_SSHF);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

double compute_mu(input in);

double normalization_condition(double mu, void *pp);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void wave_vector_grid(double *xx, input in);


// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A TWO-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx2(int xx, int yy, int x_size);


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY
// -------------------------------------------------------------------

void compute_phi(double *phi, double *xx, input in, bool verbose);

void compute_phil(double *phil, double *xx, int ll, input in);

double phixl(double yy, double xx, int ll, input in);

double phix0(double yy, double xx, input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssfHF(double *SS,  double *xx, input in);

double ssfHF(double yy, double xx, input in);


void compute_ssf(double *SS, double *SSHF,
                 double *GG, double *phi, 
		 double *xx, input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc(double *GG, double *SS, 
		  double *xx, input in);

double slfc(double yy, double xx, double SS);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

double compute_uex(double *SS, input in);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text(double *SS, double *GG, double *phi, 
		 double *xx, input in);

void write_guess(double *SS, double *GG, input in);

void read_guess(double *SS, double *GG, input in);

#endif
