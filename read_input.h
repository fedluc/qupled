#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <stdbool.h>

// -------------------------------------------------------------------
// STRUCTURE TO STORE THE INPUT PARAMETERS
// -------------------------------------------------------------------

typedef struct {

  bool guess_write;
  bool qstls_iet_static;
  char *guess_file1;
  char *guess_file2;
  char *iet_mapping;
  char *stls_guess_file;
  char *qstls_guess_file;
  char *qstls_fixed_file;
  char *qstls_iet_fixed_file;
  char *theory;
  double a_csr;
  double a_mix;
  double Theta;
  double rs;
  double dx;
  double drs;
  double dt;
  double err_min_iter;
  double mu_lo;
  double mu_hi;
  double mu;
  double xmax;
  int nl;
  int nx;
  int nrs;
  int nIter;
  int nThreads;
  
} input;

// ----------------------------------------
// FUNCTION TO READ INPUT DATA
// ----------------------------------------

void get_input(int argc, char **argv, input *in);

// -------------------------------------------------
// FUNCTION TO ASSIGN DEFAULT VALUES TO PARSER DATA
// -------------------------------------------------

void set_default_parse_opt(input *in);

// ------------------------------------------------
// FUNCTION TO COMPUTE THE NUMBER OF GRID POINTS
// ------------------------------------------------

void get_grid_size(input *in);

// ------------------------------------------------------------
// FUNCTION TO VERIFY THAT THE OPTIONS GIVEN IN INPUT ARE VALID
// ------------------------------------------------------------

void check_input(input *in);

// ------------------------------------------------
// FUNCTION TO DEBUG THE INPUT
// ------------------------------------------------

void print_input(input *in);
  
#endif
