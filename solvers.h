#ifndef SOLVERS_H
#define SOLVERS_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS EQUATIONS
// -------------------------------------------------------------------

void solve_stls(input in, bool verbose);


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE STLS-IET EQUATIONS
// -------------------------------------------------------------------

void solve_stls_iet(input in, bool verbose);


// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE QSTLS EQUATIONS
// -------------------------------------------------------------------

void solve_qstls(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTION USED TO ITERATIVELY SOLVE THE QSTLS-IET EQUATIONS
// -------------------------------------------------------------------

void solve_qstls_iet(input in, bool verbose);


#endif
