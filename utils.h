#ifndef UTILS_H
#define UTILS_H

#include "read_input.h"

// -------------------------------------------------------------------
// CONSTANTS 
// -------------------------------------------------------------------

// Tolerance to compare doubles
#define DBL_TOL 1e-14

// Maximum number of iterations for the root solvers
#define ROOT_MAX_ITER 1000

// Minimum relative error for the root solvers
#define ROOT_REL_ERR 1e-10

// Minimum absolute error for the root solvers
#define ROOT_ABS_ERR 1e-10

// Minimum relative error for the Fourier integrals
#define FOURIER_REL_ERR 1e-10

// Minimum relative error for the numerical quadratures
#define QUAD_REL_ERR 1e-5

// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A TWO-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx2(int xx, int yy, int x_size);

// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A THREE-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx3(int xx, int yy, int zz,
         int x_size, int y_size);

// -------------------------------------------------------------------
// FUNCTION USED TO GET THE SIGN OF A NUMBER
// -------------------------------------------------------------------

int get_sign(double num);

// -------------------------------------------------------------------
// FUNCTIONS USED TO READ DATA FROM TEXT FILES
// -------------------------------------------------------------------

void get_data_format_from_text(char * file_name, int *n_lines,
			       int *n_columns);

void get_data_from_text(char *file_name, int n_lines, int n_columns,
			double *data, double *xx, input *in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

double compute_internal_energy(double *SS, double *xx, input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE RADIAL DISTRIBUTION FUNCTION
// -------------------------------------------------------------------

void compute_rdf(double *gg, double *rr, double *SS,
		 double *xx, input in);

#endif
