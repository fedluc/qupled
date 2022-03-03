#ifndef QSTLS_IET_H
#define QSTLS_IET_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_qstls_iet_arrays(input in, double **psi_new);

void free_qstls_iet_arrays(double *psi_new);

// ---------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE INITIAL GUESS
// ---------------------------------------------------------------------

void initial_guess_qstls_iet(double *xx, double *SS, double *SSHF,
			     double *psi, double *phi, double *bf,
			     input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE QSTLS-IET SCHEME
// -------------------------------------------------------------------

void qstls_iet_iterations(double *SS, double *SS_new,
			  double *SSHF, double *psi,
			  double *psi_new,  double *psi_fixed_qstls,
			  double *phi, double *bf, 
			  double *xx, input in,
			  bool verbose);

void qstls_iet_update(double *SS, double *SS_new,
		      double *psi, double *psi_new,
		      input in);

// ------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FIXED COMPONENT OF THE AUXILIARY RESPONSE
// ------------------------------------------------------------------------

void init_fixed_qstls_iet_arrays(double *xx, input in, bool verbose);

void compute_adr_iet_fixed(double *xx, input in);

double adr_iet_fixed_partial_xuwl(double yy, void* pp);

double adr_iet_fixed_partial_xuw0(double yy, void* pp);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE CHANGING COMPONENT OF THE AUXILIARY RESPONSE
// ---------------------------------------------------------------------------

void compute_adr_iet(double *psi_new, double *psi, double *psi_fixed_qstls,
                     double *phi, double *SS, double *bf, double *xx,
                     input in);

double adr_iet_part1_partial(double uu, void* pp);

double adr_iet_part2_partial(double ww, void* pp);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_qstls_iet(double *SS, double *SSHF, double *psi,
			   double *phi, double *bf, double *xx, input in);



// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_qstls_iet(double *SS, double *psi, double *phi,
                          double *SSHF, double *bf, double *xx,
                          input in);

void write_text_sdr_qstls_iet(double *psi, double *phi, double *bf,
			      double *xx,  input in);

#endif
