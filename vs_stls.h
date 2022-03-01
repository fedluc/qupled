#ifndef VS_STLS_H
#define VS_STLS_H

#include "read_input.h"


// -------------------------------------------------------------------
// DATA STRUCTURES TO HANDLE MORE STATE POINTS SIMULTANEOUSLY
// -------------------------------------------------------------------

// Number of elements (field) in the vs_sp structure
const int VS_STRUCT_EL = 25;

typedef struct {

  double *rst;
  double *rstp1;
  double *rstm1;
  double *rstp2;
  double *rstm2;
  
  double *rsp1t;
  double *rsp1tp1;
  double *rsp1tm1;
  double *rsp1tp2;
  double *rsp1tm2;
  
  double *rsm1t;
  double *rsm1tp1;
  double *rsm1tm1;
  double *rsm1tp2;
  double *rsm1tm2;

  double *rsp2t;
  double *rsp2tp1;
  double *rsp2tm1;
  double *rsp2tp2;
  double *rsp2tm2;

  double *rsm2t;
  double *rsm2tp1;
  double *rsm2tm1;
  double *rsm2tp2;
  double *rsm2tm2;
 
} vs_struct;

// -------------------------------------------------------------------
// DATA STRUCTURE USED TO ALLOCATE THE THERMO VS-STLS ARRAYS
// -------------------------------------------------------------------

typedef struct {

  double *rsu;
  double *rsa;
  
} vs_stls_pointers;


// -------------------------------------------------------------------
// FUNCTION USED TO LOOP OVER THE VS_STRUCT DATA STRUCTURES
// -------------------------------------------------------------------

double *get_el(vs_struct sp, int ii);

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_vs_stls_struct_arrays(input in, stls_pointers pp, vs_struct *xx,
				vs_struct *phi, vs_struct *SS, vs_struct *SSHF,
				vs_struct *GG, vs_struct *GG_new, int el);

void alloc_vs_stls_thermo_arrays(input in, vs_stls_pointers pp,
				 vs_struct *rsu, vs_struct *rsa, int el);
  
void alloc_vs_stls_thermo_arrays_elem(input in, vs_stls_pointers *pp);
  
void alloc_vs_stls_one_array(double *pp, vs_struct *vs_arr, int el);

void free_vs_stls_struct_arrays(vs_struct xx, vs_struct phi, vs_struct SS,
				vs_struct SSHF, vs_struct GG, vs_struct GG_new);

void free_vs_stls_thermo_arrays(vs_struct rsu, vs_struct rsa);

void free_vs_stls_one_array(vs_struct vs_arr);

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_vs_stls_arrays(input *in, input *vs_in, vs_struct xx,
			       vs_struct rsa, bool verbose);

void init_tmp_vs_stls_arrays(input *vs_in, vs_struct xx,
			     vs_struct phi, vs_struct SSHF,
			     bool verbose);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rsp, input *in);

// ---------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE INITIAL GUESS
// ---------------------------------------------------------------------

void initial_guess_vs_stls(vs_struct xx, vs_struct SS, vs_struct SSHF,
			   vs_struct GG, vs_struct GG_new, vs_struct phi,
			   input *vs_in);

// ---------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE VS-STLS SCHEME
// ---------------------------------------------------------------------

double vs_stls_thermo_iterations(vs_struct xx, vs_struct rsu, vs_struct rsa,
				 input *vs_in, bool verbose);

double vs_stls_thermo_err(double alpha, input *vs_in);

void vs_stls_thermo_update(double alpha, input *vs_in);

void vs_stls_struct_iterations(vs_struct SS, vs_struct SSHF,
			       vs_struct GG, vs_struct GG_new,
			       vs_struct phi, vs_struct xx,
			       input *vs_in, bool verbose);

double vs_stls_struct_err(vs_struct GG, vs_struct GG_new, input *vs_in);

void vs_stls_struct_update(vs_struct GG, vs_struct GG_new, input *vs_in);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

void compute_chemical_potential_vs_stls(input *vs_in);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_vs_stls(vs_struct SS, vs_struct SSHF, vs_struct GG, vs_struct phi,
		    vs_struct xx, input *vs_in);

void compute_ssf_HF_vs_stls(vs_struct SS, vs_struct xx, input *vs_in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_vs_stls(vs_struct GG, vs_struct SS,
			  vs_struct xx, input *vs_in);

// -----------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY RESPONSE
// -----------------------------------------------------------------------

void compute_idr_vs_stls(vs_struct phi, vs_struct xx, input *vs_in);
  
// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE PARAMETER FOR THE CSR RULE
// -------------------------------------------------------------------

double compute_alpha(vs_struct xx, vs_struct rsu, vs_struct rsa, input *vs_in);


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE INTEGRAND FOR THE FREE ENERGY
// -------------------------------------------------------------------

void compute_rsu(vs_struct xx, vs_struct rsu, vs_struct rsa,
		 input *vs_in, bool verbose);

void compute_rsu_blocks(vs_struct xx, vs_struct rsu, vs_struct rsa,
			input *vs_in, int start, int end,
			int step, bool verbose);
    
// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_vs_stls(double *rsu, double *rsp, input in);

void read_guess_vs_stls(vs_struct SS, vs_struct GG, input *vs_in);
  
#endif

