#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "stls.h"
#include "stls_hnc.h"
#include "qstls.h"

// ----------------------------------------
// COMMAND LINE PARSER
// ----------------------------------------

// Documentation
static char doc[] =
  "The documentation is available at https://github.com/fedluc/STLS";



// Optional arguments
static struct argp_option options[] = {
  {"Theta",    't', "1.0", 0,
   "Quantum degeneracy parameter"}, 
  {"rs",    'r', "1.0", 0,
   "Quantum coupling parameter"},
  {"xcut",  'x', "20.48", 0,
   "Cutoff for wave-vector grid"},
  {"dx",  'd', "0.01", 0,
   "Resolution for wave-vector grid"},
  {"nl",  'l', "128", 0,
   "Number of Matsubara frequencies"},
  {"iter",  'i', "1000", 0,
   "Maximum number of iterations"},
  {"errIter", 'e', "1e-5", 0,
   "Minimum error for convergence in the iterations"},
  {"mix", 'm', "0.1", 0,
   "Mixing parameter for iterative solution"},
  {"mg", 'g', "-10,10", 0,
   "Initial guess for chemical potential"},
  {"phi", 'p', "NO_FILE", 0,
   "Load density response from PHI_FILE"},
  {"uex", 'u', "NO_FILE", 0,
   "Compute internal energy from data in SSF_FILE"},
  {"sol", 's', "STLS",0,
   "Theory to be solved"},
  { 0 }
};

// Structure to communicate between main and parser
struct arguments
{

  char *phi_file;
  char *ssf_file;
  char *theory; 
  float Theta;
  float rs;
  float dx;
  float err_min_iter;
  float err_min_int;
  float a_mix;
  float mu_lo;
  float mu_hi;
  float xmax;
  int nl;
  int nIter;
  
};


// Single option parser
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{

  char *value;
  struct arguments *arguments = state->input;
  
  switch (key)
    {

    case 'd':
      arguments->dx = atof(arg);
      break;
    case 'e':
      arguments->err_min_iter = atof(arg);
      break;
    case 'g':
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
    case 'i':
      arguments->nIter = atoi(arg);
      break;
    case 'l':
      arguments->nl = atoi(arg);
      break;
    case 'm':
      arguments->a_mix = atof(arg);
      break;
    case 'p':
      arguments->phi_file = arg;
      break;
    case 'r':
      arguments->rs = atof(arg);
      break;
    case 't':
      arguments->Theta = atof(arg);
      break;
    case 's':
      arguments->theory = arg;
      break;
    case 'u':
      arguments->ssf_file = arg;
      break;
    case 'x':
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
  arguments.phi_file = "NO_FILE"; // File with density response
  arguments.ssf_file = "NO_FILE"; // File with static structure factor
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

  // Parse command line
  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  // Fill input structure
  input in;
  in.phi_file = arguments.phi_file;
  in.ssf_file = arguments.phi_file;
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

  // Solve theory specified in input
  clock_t tic = clock();
  if (strcmp(arguments.theory, "STLS") == 0)
    solve_stls(in, true, NULL, NULL, NULL, NULL, NULL, NULL);
  else if (strcmp(arguments.theory, "STLS-HNC") == 0)
    solve_stls_hnc(in, true);
  else if (strcmp(arguments.theory, "QSTLS") == 0)
    solve_qstls(in, true);
  else
    printf("Error: unknown theory to be solved. Choose between: STLS, STLS-HNC and QSTLS\n");
  clock_t toc = clock();
  printf("Solution complete. Elapsed time: %f seconds\n", ((double)toc - (double)tic) / CLOCKS_PER_SEC);

  return 0;

}

