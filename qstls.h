#ifndef QSTLS_H
#define QSTLS_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------

void alloc_qstls_arrays(input in, double **psi, double **psi_fixed);

void free_qstls_arrays(double *psi, double *psi_fixed);

// -------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE INITIAL GUESS
// -------------------------------------------------------------------

void initial_guess_qstls(double *xx, double *SS, double *SSHF,
			 double *psi, double *phi, input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE STLS SCHEME
// -------------------------------------------------------------------

void qstls_iterations(double *SS, double *SS_new,
		      double *SSHF, double *psi,
		      double *psi_fixed, double *phi,
		      double *xx, input in,
		      bool verbose);
  
// ------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FIXED COMPONENT OF THE AUXILIARY RESPONSE
// ------------------------------------------------------------------------

void init_fixed_qstls_arrays(double *psi_fixed, double *xx,
			     input in, bool verbose);

void compute_adr_fixed(double *psi_fixed, double *xx, input in);

double adr_fixed_lev1_partial_xwl(double qq, void* pp);

double adr_fixed_lev1_partial_xw0(double qq, void* pp);

double adr_fixed_lev2_partial_xwq0(double tt, void* pp);

double adr_fixed_lev2_partial_xwql(double tt, void* pp);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE CHANGING COMPONENT OF THE AUXILIARY RESPONSE
// ---------------------------------------------------------------------------

void compute_adr(double *psi, double *psi_fixed, double *SS,
		 double *xx, input in);

double adr_partial(double ww, void* pp);


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_qstls(double *SS, double *SSHF, double *psi,
		       double *phi, double *xx, input in);


// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_qstls(double *SS, double *psi, double *phi, 
			double *SSHF, double *xx, input in);

void write_text_slfc_qstls(double *psi, double *phi,
			   double *xx,  input in);

void write_text_sdr_qstls(double *psi, double *phi,
			  double *xx,  input in);

void write_text_adr(double *psi, input in);

void write_guess_qstls(double *SS, double *psi, input in);

void read_guess_qstls(double *SS, double *psi, input in);

void write_fixed_qstls(double *psi_fixed, input in);

void read_fixed_qstls(double *psi_fixed, input in);

void check_guess_qstls(int nx, double dx, double xmax,
		       int nl, double Theta, input in,
		       size_t it_read, size_t it_expected,
		       FILE *fid, bool check_grid,
		       bool check_items, bool check_eof);

#endif
