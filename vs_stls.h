
#ifndef VS_STLS_H
#define VS_STLS_H

#include "read_input.h"

// ---------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE VS-STLS SCHEME
// ---------------------------------------------------------------------

double vs_stls_thermo_iterations(double *xx, double *uint,
				 double *rsArray, input in,
				 bool verbose);

void fill_internal_energy_array(double *xx, double *uint,
				double *rsArray, input in,
				bool verbose);

void vs_stls_struct_iterations(double *SS, double *SSHF,
			       double *GG, double *GG_new,
			       double *phi, double *xx,
			       input in, bool verbose);
  
// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_vs_stls_arrays(input in, double **uint, double **rsArray);

void free_vs_stls_arrays(double *uint, double *rsArray);

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_vs_stls_arrays(input *in, double *xx, double *rsArray, bool verbose);


void init_state_point_vs_stls_arrays(input *in, double *xx,
				     double *phi, double *SSHF,
				     bool verbose);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rsArray, input *in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_vs_slfc(double *GG, double *SS, double *xx, input in);


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE PARAMETER FOR THE CSR RULE
// -------------------------------------------------------------------

double get_alpha(double *xx, double *uint, double *rsArray, input in);


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE INTERNAL ENERGY FOR THERMODYNAMIC
// INTEGRATION
// -------------------------------------------------------------------

void fill_internal_energy_array(double *xx, double *uint, double *rsArray,
				input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FREE ENERGY
// -------------------------------------------------------------------

double compute_free_energy(double *uint, double *rsArray, input in);

double fex(double rs, void* pp);


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_vs_stls(double *uint, double *rsArray, input in);
  
#endif

