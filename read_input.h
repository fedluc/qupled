#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <stdbool.h>

// -------------------------------------------------------------------
// CONSTANTS 
// -------------------------------------------------------------------

// Default string for when no file name is passed in input
#define NO_FILE_STR "\0"

// -------------------------------------------------------------------
// STRUCTURE TO STORE THE INPUT PARAMETERS
// -------------------------------------------------------------------

typedef struct {

  bool qstls_iet_static;
  bool vs_solve_csr;
  char *dyn_adr_file;
  char *guess_file1;
  char *guess_file2;
  char *mode;
  char *iet_mapping;
  char *stls_guess_file;
  char *qstls_guess_file;
  char *qstls_fixed_file;
  char *qstls_iet_fixed_file;
  char *theory;
  char *vs_thermo_file;
  double a_mix;
  double Theta;
  double rs;
  double dx;
  double dyn_dW;
  double dyn_Wmax;
  double dyn_xtarget;
  double err_min_iter;
  double mu_lo;
  double mu_hi;
  double mu;
  double vs_alpha;
  double vs_a_mix;
  double vs_drs;
  double vs_dt;
  double vs_err_min_iter;
  double xmax;
  int nl;
  int nx;
  int nIter;
  int nThreads;
  int nW;
  int vs_nrs;
  
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
