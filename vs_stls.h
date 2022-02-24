#ifndef VS_STLS_H
#define VS_STLS_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_vs_stls_arrays(input in, int nrs, double **uint, double **rsArray);

void free_vs_stls_arrays(double *uint, double *rsArray);

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_vs_stls_arrays(input *in, int nrs, double *xx, double *rsArray, bool verbose);


void init_state_point_vs_stls_arrays(input *in, double *xx,
				     double *phi, double *SSHF,
				     bool verbose);

void vs_stls_iterations(double *SS, double *SSHF,
			double *GG, double *GG_new,
			double *phi, double *xx,
			double drs, double alpha,
			input in, bool verbose);

void compute_vs_slfc(double *GG, double *SS, double *xx,
		     double drs, double alpha, input in);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rsArray, int nrs, input *in);



// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FREE ENERGY
// -------------------------------------------------------------------

double compute_free_energy(double *uint, double *rsArray, input in, int nrs);

double fex(double rs, void* pp);

#endif

