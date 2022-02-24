#ifndef UTILS_H
#define UTILS_H

#include "read_input.h"

// -------------------------------------------------------------------
// CONSTANTS FOR ROOT SOLVERS, MINIMIZATIONS AND QUADRATURES
// -------------------------------------------------------------------

// Maximum number of iterations for the root solvers
static const int ROOT_MAX_ITER = 1000;

// Minimum relative error for the root solvers
static const double ROOT_REL_ERR = 1e-10;

// Minimum absolute error for the root solvers
static const double ROOT_ABS_ERR = 1e-10;

// Minimum relative error for the Fourier integrals
static const double FOURIER_REL_ERR = 1e-10;

// Minimum relative error for the numerical quadratures
static const double QUAD_REL_ERR = 1e-5;

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
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

double compute_internal_energy(double *SS, double *xx, input in);

double uex(double yy, void *pp);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FREE ENERGY
// -------------------------------------------------------------------

double compute_free_energy(double *rsu, double *rsp, input in);

double fex(double rs, void* pp);


#endif
