#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <argp.h>
#include <math.h>
#include <omp.h>
#include "solvers.h"

// ----------------------------------------
// COMMAND LINE PARSER
// ----------------------------------------

// Documentation
static char doc[] =
  "The documentation is available at https://github.com/fedluc/STLS";

// Non-ascii characters to avoid using short options in the parser
#define ARGUMENT_THETA_SHORT 0x80
#define ARGUMENT_RS_SHORT 0x81
#define ARGUMENT_XMAX_SHORT 0x82
#define ARGUMENT_DX_SHORT 0x83
#define ARGUMENT_NL_SHORT 0x84
#define ARGUMENT_ITER_SHORT 0x85
#define ARGUMENT_MIN_ERR_SHORT 0x86
#define ARGUMENT_MIX_SHORT 0x87
#define ARGUMENT_MU_GUESS_SHORT 0x88
#define ARGUMENT_STLS_GUESS_SHORT 0x89
#define ARGUMENT_THEORY_SHORT 0x90
#define ARGUMENT_OMP_SHORT 0x91
#define ARGUMENT_QSTLS_GUESS_SHORT 0x92
#define ARGUMENT_QSTLS_FIXED_SHORT 0x93
#define ARGUMENT_QSTLS_IET_FIXED_SHORT 0x94
#define ARGUMENT_QSTLS_IET_STATIC_SHORT 0x95


// Optional arguments
static struct argp_option options[] = {
  {"Theta", ARGUMENT_THETA_SHORT, "1.0", 0,
   "Quantum degeneracy parameter"},
  {"rs", ARGUMENT_RS_SHORT, "1.0", 0,
   "Quantum coupling parameter"},
  {"xmax", ARGUMENT_XMAX_SHORT, "20.48", 0,
   "Cutoff for wave-vector grid"},
  {"dx", ARGUMENT_DX_SHORT, "0.01", 0,
   "Resolution for wave-vector grid"},
  {"nl", ARGUMENT_NL_SHORT, "128", 0,
   "Number of Matsubara frequencies"},
  {"iter", ARGUMENT_ITER_SHORT, "1000", 0,
   "Maximum number of iterations"},
  {"min-err", ARGUMENT_MIN_ERR_SHORT, "1e-5", 0,
   "Minimum error for convergence in the iterations"},
  {"mix", ARGUMENT_MIX_SHORT, "0.1", 0,
   "Mixing parameter for iterative solution"},
  {"mu-guess", ARGUMENT_MU_GUESS_SHORT, "-10,10", 0,
   "Initial guess for chemical potential"},
  {"stls-guess", ARGUMENT_STLS_GUESS_SHORT, "NO_FILE", 0,
   "Load initial guess from file for the stls and stls-iet schemes"},
  {"qstls-guess", ARGUMENT_QSTLS_GUESS_SHORT, "NO_FILE", 0,
   "Load initial guess from file for the qstls and qstls-iet schemes"},
  {"qstls-fix", ARGUMENT_QSTLS_FIXED_SHORT, "NO_FILE", 0,
   "Load fixed component of the density response function from file for the qslts scheme"},
  {"qstls-iet-fix", ARGUMENT_QSTLS_IET_FIXED_SHORT, "NO_FILE", 0,
   "Load fixed component of the density response function from file for the qslts-iet scheme"},
  {"qstls-iet-static", ARGUMENT_QSTLS_IET_STATIC_SHORT, "0", 0,
   "Use static approximation to compute the auxilliary density response for the qstls-iet scheme"},
  {"theory", ARGUMENT_THEORY_SHORT, "STLS",0,
   "Scheme to be solved"},
  {"omp", ARGUMENT_OMP_SHORT, "1",0,
   "Number of omp threads to use in the solution"},
  { 0 }
};

// Structure to communicate between main and parser
struct arguments
{

  char *stls_guess_file;
  char *qstls_guess_file;
  char *qstls_fixed_file;
  char *qstls_iet_fixed_file;
  char *theory; 
  double Theta;
  double rs;
  double dx;
  double err_min_iter;
  double err_min_int;
  double a_mix;
  double mu_lo;
  double mu_hi;
  double xmax;
  int nl;
  int nIter;
  int nThreads;
  int qstls_iet_static;

};


// Single option parser
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{

  char *value;
  struct arguments *arguments = state->input;
  
  switch (key)
    {

    case ARGUMENT_DX_SHORT:
      arguments->dx = atof(arg);
      break;
    case  ARGUMENT_MIN_ERR_SHORT:
      arguments->err_min_iter = atof(arg);
      break;
    case  ARGUMENT_STLS_GUESS_SHORT:
      arguments->stls_guess_file = arg;
      break;
    case  ARGUMENT_QSTLS_GUESS_SHORT:
      arguments->qstls_guess_file = arg;
      break;
    case  ARGUMENT_QSTLS_FIXED_SHORT:
      arguments->qstls_fixed_file = arg;
      break;
    case  ARGUMENT_QSTLS_IET_FIXED_SHORT:
      arguments->qstls_iet_fixed_file = arg;
      break;
    case  ARGUMENT_QSTLS_IET_STATIC_SHORT:
      arguments->qstls_iet_static = atoi(arg);
      break;
    case  ARGUMENT_MU_GUESS_SHORT:
      value = strtok(NULL, ",");
      if(value != NULL ) {
	arguments->mu_lo = atof(value);
      }
      else exit(EXIT_FAILURE);
      value = strtok(NULL, ",");
      if(value != NULL ) {
	arguments->mu_hi = atof(value);
      }
      else exit(EXIT_FAILURE);
      break;
    case ARGUMENT_ITER_SHORT:
      arguments->nIter = atoi(arg);
      break;
    case ARGUMENT_NL_SHORT:
      arguments->nl = atoi(arg);
      break;
    case ARGUMENT_MIX_SHORT:
      arguments->a_mix = atof(arg);
      break;
    case ARGUMENT_OMP_SHORT:
      arguments->nThreads = atoi(arg);
      break;
    case ARGUMENT_RS_SHORT:
      arguments->rs = atof(arg);
      break;
    case ARGUMENT_THETA_SHORT:
      arguments->Theta = atof(arg);
      break;
    case ARGUMENT_THEORY_SHORT:
      arguments->theory = arg;
      break;
    case  ARGUMENT_XMAX_SHORT:
      arguments->xmax = atof(arg);
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num > 0) 
        argp_usage (state);

      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, 0, doc };

// ----------------------------------------
// MAIN
// ----------------------------------------

int main (int argc, char **argv){

  struct arguments arguments;

  // Default values for optional arguments
  arguments.stls_guess_file = "NO_FILE"; // File with initial guess for STLS and STLS-IET schemes
  arguments.qstls_guess_file = "NO_FILE"; // File with initial guess for QSTLS and QSTLS-IET schemes
  arguments.qstls_fixed_file = "NO_FILE"; // File with fixed component of the density response for the QSTLS scheme
  arguments.qstls_iet_fixed_file = "NO_FILE"; // File with fixed component of the density response for the QSTLS-IET scheme
  arguments.qstls_iet_static = 0; // Static approximation for the auxiliary density response for the QSTLS-IET scheme (0 = off, 1 = on)
  arguments.Theta = 1.0; // Quantum degeneracy parameter
  arguments.rs = 1.0; // Quantum coupling parameter
  arguments.dx = 0.01; // Wave-vector grid resolution
  arguments.err_min_iter = 1e-5; // Minimum error for convergence in the iterative procedure
  arguments.a_mix = 0.1; // Mixing parameter for iterative procedure
  arguments.mu_lo = -10; // Initial guess for chemical potential (low bound)
  arguments.mu_hi = 10; // Initial guess for chemical potential (high bound)
  arguments.xmax = 20.48; // Cutoff for wave-vector grid
  arguments.nl = 128; // Number of Matsubara frequencies
  arguments.nIter = 1000; // Number of iterations
  arguments.theory = "STLS"; // Theory to solve
  arguments.nThreads = 1; // Number of OMP threads to use in the solution

  // Parse command line
  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  // Fill input structure
  input in;
  in.stls_guess_file = arguments.stls_guess_file;
  in.qstls_guess_file = arguments.qstls_guess_file;
  in.qstls_fixed_file = arguments.qstls_fixed_file;
  in.qstls_iet_fixed_file = arguments.qstls_iet_fixed_file;
  in.qstls_iet_static = arguments.qstls_iet_static;
  in.Theta = arguments.Theta;
  in.rs = arguments.rs;
  in.dx = arguments.dx; 
  in.err_min_iter = arguments.err_min_iter;
  in.a_mix = arguments.a_mix;
  in.mu_lo = arguments.mu_lo;
  in.mu_hi = arguments.mu_hi;
  in.xmax = arguments.xmax;
  in.nl = arguments.nl;
  in.nIter = arguments.nIter;
  in.nx = (int)floor(in.xmax/in.dx);
  in.theory = arguments.theory;
  if (strcmp(arguments.theory, "STLS") == 0) in.theory_id = 1;
  else if (strcmp(arguments.theory, "STLS-IET-HNC") == 0) in.theory_id = 2;
  else if (strcmp(arguments.theory, "STLS-IET-IOI") == 0) in.theory_id = 3;
  else if (strcmp(arguments.theory, "STLS-IET-LCT") == 0) in.theory_id = 4;
  else if (strcmp(arguments.theory, "STLS-RIET-LCT") == 0) in.theory_id = 5;
  else if (strcmp(arguments.theory, "QSTLS") == 0) in.theory_id = 6;
  else if (strcmp(arguments.theory, "QSTLS-IET-HNC") == 0) in.theory_id = 7;
  else if (strcmp(arguments.theory, "QSTLS-IET-IOI") == 0) in.theory_id = 8;
  else if (strcmp(arguments.theory, "QSTLS-IET-LCT") == 0) in.theory_id = 9;
  else if (strcmp(arguments.theory, "QSTLS-RIET-LCT") == 0) in.theory_id = 10;
  else {
    printf("Error: unknown theory to be solved." 
	   "Choose between: STLS, STLS-IET-HNC, STLS-IET-IOI,"
	   "STLS-IET-LCT, STLS-RIET-LCT, QSTLS, QSTLS-IET-HNC,"
	   "QSTLS-IET-IOI, QSTLS-IET-LCT and QSTLS-RIET-LCT\n");
    exit(EXIT_FAILURE);
  }

  // Set number of threads for parallel calculations
  omp_set_num_threads(arguments.nThreads);


  // Start timing
  double tic = omp_get_wtime();

  // Solve theory specified in input
  if (in.theory_id == 1)
    solve_stls(in, true);
  else if (in.theory_id == 2 || 
	   in.theory_id == 3 ||
	   in.theory_id == 4 ||
           in.theory_id == 5)
    solve_stls_iet(in, true);
  else if (in.theory_id == 6)
    solve_qstls(in, true);
  else if (in.theory_id == 7 ||
	   in.theory_id == 8 ||
	   in.theory_id == 9 ||
	   in.theory_id == 10)
    solve_qstls_iet(in, true);

  // End timing 
  double toc = omp_get_wtime();

  // Print conclusion message on screen
  printf("Solution complete. Elapsed time: %f seconds\n", toc - tic);



  return 0;

}
