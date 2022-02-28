#ifndef VS_STLS_H
#define VS_STLS_H

#include "read_input.h"


// -------------------------------------------------------------------
// DATA STRUCTURES TO HANDLE MORE STATE POINTS SIMULTANEOUSLY
// -------------------------------------------------------------------

// Number of elements (field) in the vs_sp structure
const int VS_SP_EL = 5;

typedef struct {

  double *in;
  double *rsp1;
  double *rsm1;
  double *rsp2;
  double *rsm2;
  double *tp1;
  double *tm1;
  double *tp2;
  double *tm2;
  
} vs_sp;

// -------------------------------------------------------------------
// DATA STRUCTURE USED TO ALLOCATE THE THERMO VS-STLS ARRAYS
// -------------------------------------------------------------------

typedef struct {

  double *rsu;
  double *rsa;
  
} vs_stls_pointers;


// -------------------------------------------------------------------
// FUNCTION USED TO LOOP OVER THE VS_SP DATA STRUCTURES
// -------------------------------------------------------------------

double *get_el(vs_sp sp, int ii);

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_vs_stls_struct_arrays(input in, stls_pointers pp, vs_sp *xx,
				vs_sp *phi, vs_sp *SS, vs_sp *SSHF,
				vs_sp *GG, vs_sp *GG_new, int el);

void alloc_vs_stls_thermo_arrays(input in, vs_stls_pointers pp,
				 vs_sp *rsu, vs_sp *rsa, int el);
  
void alloc_vs_stls_thermo_arrays_elem(input in, vs_stls_pointers *pp);
  
void alloc_vs_stls_one_array(double *pp, vs_sp *vs_arr, int el);

void free_vs_stls_struct_arrays(vs_sp xx, vs_sp phi, vs_sp SS,
				vs_sp SSHF, vs_sp GG, vs_sp GG_new);

void free_vs_stls_thermo_arrays(vs_sp rsu, vs_sp rsa);

void free_vs_stls_one_array(vs_sp vs_arr);

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_vs_stls_arrays(input *in, input *vs_in, vs_sp xx,
			       vs_sp rsa, bool verbose);

void init_state_point_vs_stls_arrays(input *vs_in, vs_sp xx,
				     vs_sp phi, vs_sp SSHF,
				     bool verbose);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rsp, input *in);

// ---------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE INITIAL GUESS
// ---------------------------------------------------------------------

void initial_guess_vs_stls(vs_sp xx, vs_sp SS, vs_sp SSHF,
			   vs_sp GG, vs_sp GG_new, vs_sp phi,
			   input *vs_in);

// ---------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE VS-STLS SCHEME
// ---------------------------------------------------------------------

double vs_stls_thermo_iterations(vs_sp xx, vs_sp rsu, vs_sp rsa,
				 input *vs_in, bool verbose);

double vs_stls_thermo_err(double alpha, input *vs_in);

void vs_stls_thermo_update(double alpha, input *vs_in);

void vs_stls_struct_iterations(vs_sp SS, vs_sp SSHF,
			       vs_sp GG, vs_sp GG_new,
			       vs_sp phi, vs_sp xx,
			       input *vs_in, bool verbose);

double vs_stls_struct_err(vs_sp GG, vs_sp GG_new, input *vs_in);

void vs_stls_struct_update(vs_sp GG, vs_sp GG_new, input *vs_in);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

void compute_chemical_potential_vs_stls(input *vs_in);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_vs_stls(vs_sp SS, vs_sp SSHF, vs_sp GG, vs_sp phi,
		    vs_sp xx, input *vs_in);

void compute_ssf_HF_vs_stls(vs_sp SS, vs_sp xx, input *vs_in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_vs_stls(vs_sp GG, vs_sp SS,
			  vs_sp xx, input *vs_in);

// -----------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY RESPONSE
// -----------------------------------------------------------------------

void compute_idr_vs_stls(vs_sp phi, vs_sp xx, input *vs_in);
  
// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE PARAMETER FOR THE CSR RULE
// -------------------------------------------------------------------

double compute_alpha(vs_sp xx, vs_sp rsu, vs_sp rsa, input *vs_in);


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE INTEGRAND FOR THE FREE ENERGY
// -------------------------------------------------------------------

void compute_rsu(vs_sp xx, vs_sp rsu, vs_sp rsa,
		 input *vs_in, bool verbose);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_vs_stls(double *rsu, double *rsp, input in);

void read_guess_vs_stls(vs_sp SS, vs_sp GG, input *vs_in);
  
#endif

