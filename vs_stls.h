#ifndef VS_STLS_H
#define VS_STLS_H

#include "read_input.h"


// -------------------------------------------------------------------
// STRUCTURE TO STORE THE RESULTS FOR MULTIPLE STATE POINTS
// -------------------------------------------------------------------

// Number of elements (field) in the vs_sp structure
const int VS_SP_EL = 5;

typedef struct {

  double *in;
  double *rsp1;
  double *rsm1;
  double *rsp2;
  double *rsm2;
  
} vs_sp;

typedef struct {

  double *xx;
  double *phi;
  double *SS;
  double *SSHF;
  double *GG;
  double *GG_new;
  
} stls_pointers;

// ---------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE VS-STLS SCHEME
// ---------------------------------------------------------------------

double vs_stls_thermo_iterations(vs_sp xx, double *rsu,
				 double *rsp, input *vs_in,
				 bool verbose);

void vs_stls_struct_iterations(vs_sp SS, vs_sp SSHF,
			       vs_sp GG, vs_sp GG_new,
			       vs_sp phi, vs_sp xx,
			       input *vs_in, bool verbose);
  
// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_vs_stls_arrays(input in, stls_pointers pp, vs_sp *xx,
			  vs_sp *phi, vs_sp *SS, vs_sp *SSHF,
			  vs_sp *GG, vs_sp *GG_new, double **rsu,
			  double **rs_arr, int el);

void free_vs_stls_arrays_new(vs_sp xx, vs_sp phi, vs_sp SS, vs_sp SSHF,
			     vs_sp GG, vs_sp GG_new, double *rsu,
			     double *rs_arr);

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_vs_stls_arrays(input *in, input *vs_in, vs_sp xx,
			       double *rs_arr, bool verbose);

void init_state_point_vs_stls_arrays(input *vs_in, vs_sp xx,
				     vs_sp phi, vs_sp SSHF,
				     bool verbose);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rsp, input *in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_vs_slfc(vs_sp GG, vs_sp SS,
		     vs_sp xx, input *vs_in);


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE PARAMETER FOR THE CSR RULE
// -------------------------------------------------------------------

double compute_alpha(vs_sp xx, double *rsu,
		     double *rsp, input *vs_in);


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE INTEGRAND FOR THE FREE ENERGY
// -------------------------------------------------------------------

void compute_rsu(vs_sp xx, double *rsu, double *rsp,
		 input *vs_in, bool verbose);


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_vs_stls(double *rsu, double *rsp, input in);
  
#endif

