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

void init_fixed_qstls_arrays(double *psi_fixed, double *xx,
			     input in, bool verbose);

// ---------------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE AUXILIARY DENSITY RESPONSE
// ---------------------------------------------------------------------------

void compute_adr(double *psi, double *psi_fixed, double *SS,
		 double *xx, input in);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_slfc_qstls(double *psi, double *phi,
			   double *xx,  input in);
void write_text_adr(double *psi, input in);

void write_restart_qstls(double *SS, double *psi, input in);

void read_restart_qstls(double *SS, double *psi, input in);

void check_restart_qstls(int nx, double dx, double xmax,
			 int nl, double Theta, input in,
			 size_t it_read, size_t it_expected,
			 FILE *fid, bool check_grid,
			 bool check_items, bool check_eof);

#endif
