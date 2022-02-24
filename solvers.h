#ifndef SOLVERS_H
#define SOLVERS_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTION USED TO SOLVE THE STLS SCHEME
// -------------------------------------------------------------------

void solve_stls(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTION USED TO SOLVE THE VS-STLS SCHEME
// -------------------------------------------------------------------

void solve_vs_stls(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTION USED TO SOLVE THE STLS-IET SCHEME
// -------------------------------------------------------------------

void solve_stls_iet(input in, bool verbose);


// -------------------------------------------------------------------
// FUNCTION USED TO SOLVE THE QSTLS SCHEME
// -------------------------------------------------------------------

void solve_qstls(input in, bool verbose);

// -------------------------------------------------------------------
// FUNCTION USED TO SOLVE THE QSTLS-IET SCHEME
// -------------------------------------------------------------------

void solve_qstls_iet(input in, bool verbose);


#endif
