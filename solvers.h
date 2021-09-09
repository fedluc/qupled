#ifndef SOLVERS_H
#define SOLVERS_H

#include <stdbool.h>

// -------------------------------------------------------------------
// STRUCTURE TO STORE THE INPUT PARAMETERS
// -------------------------------------------------------------------

typedef struct {

  char *guess_file;
  char *theory;
  double Theta;
  double rs;
  double dx;
  double err_min_iter;
  double a_mix;
  double mu_lo;
  double mu_hi;
  double mu;
  double xmax;
  int nl;
  int nx;
  int nIter;
  int theory_id;


} input;

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
