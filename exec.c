#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include <omp.h>
#include "stls.h"

// ----------------------------------------
// COMMAND LINE PARSER
// ----------------------------------------

// Documentation
static char doc[] =
  "stls solves the classical STLS approach as defined in S. Tanaka\n"
  "and S. Ichimaru, J. Phys. Soc. Jpn. 55, 2278 (1986). The state \n"
  "point of interest is defined via the quantum degeneracy parameter\n"
  "and via the quantum coupling parameters. The equation for the \n"
  "chemical potential is solved via bisection method for which two \n"
  "initial guesses must be provided in input via the option -g.\n" 
  "The STLS approach is solved iteratively on a wave-vector grid \n"
  "extending from 0 to agiven cutoff. The grid resolution and cutoff\n"
  "are specified in input together with number of Matsubara frequencies\n"
  " necessary for the calculation of the static structure factor. The \n"
  "iterative solution employs mixing [K.-C. Ng, J. Chem. Phys. 61 2680 \n"
  "(1974)] and is assumed to have converged once the condition \n"
  "||G_{i}(x) - G_{i-1}(x)|| < epsilon is satisfied. Here G(x) is the \n"
  "static local field correction and epsilon is a tolerance specified in\n" 
  "input. The output of the code consists of:\n"
  "    - One text file with the static structure factor\n"
  "    - One text file with the static local field correction\n"
  "    - One binary file with the density response. Since the density\n"
  "      response depends only on Theta, this file can be stored and\n"
  "      provided in input for subsequent solutions of the STLS approach\n"
  "      with the same Theta (see option -p or --phi). It should be noted\n"
  "      that, if the option -p is used, the values of the quantum\n"
  "      degeneracy parameter, of the grid resolution and of the grid\n"
  "      cutoff specified in input will be overwritten by the values\n" 
  "      contained in the density response file provided in input";



// Optional arguments
static struct argp_option options[] = {
  {"Theta",    't', "1.0", 0,
   "Quantum degeneracy parameter"}, 
  {"rs",    'r', "1.0", 0,
   "Quantum coupling parameter"},
  {"xcut",  'x', "50", 0,
   "Cutoff for wave-vector grid"},
  {"dx",  'd', "0.01", 0,
   "Resolution for wave-vector grid"},
  {"nl",  'l', "1000", 0,
   "Number of Matsubara frequencies"},
  {"iter",  'i', "1000", 0,
   "Maximum number of iterations"},
  {"err", 'e', "1e-5", 0,
   "Minimum error for convergence"},
  {"mix", 'm', "0.1", 0,
   "Mixing parameter for iterative solution"},
  {"mg", 'g', "-10,10", 0,
   "Initial guess for chemical potential"},
  {"phi", 'p', "NO_FILE", 0,
   "Load density response from PHI_FILE"},
  {"uex", 'u', "NO_FILE", 0,
   "Compute internal energy from data in SSF_FILE"},
  { 0 }
};

// Structure to communicate between main and parser
struct arguments
{

  char *phi_file;
  char* ssf_file;
  double Theta;
  double rs;
  double dx;
  double err_min;
  double a_mix;
  double mu_lo;
  double mu_hi;
  double xmax;
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
      arguments->err_min = atof(arg);
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
  arguments.err_min = 1e-5; // Minimum error for convergence
  arguments.a_mix = 0.1; // Mixing parameter for iterative procedure
  arguments.mu_lo = -10; // Initial guess for chemical potential (low bound)
  arguments.mu_hi = 10; // Initial guess for chemical potential (high bound)
  arguments.xmax = 50; // Cutoff for wave-vector grid
  arguments.nl = 1000; // Number of Matsubara frequencies
  arguments.nIter = 1000; // Number of iterations 

  // Parse command line
  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  // Fill input structure
  input in;
  in.phi_file = arguments.phi_file;
  in.ssf_file = arguments.phi_file;
  in.Theta = arguments.Theta;
  in.rs = arguments.rs;
  in.dx = arguments.dx;
  in.err_min = arguments.err_min;
  in.a_mix = arguments.a_mix;
  in.mu_lo = arguments.mu_lo;
  in.mu_hi = arguments.mu_hi;
  in.xmax = arguments.xmax;
  in.nl = arguments.nl;
  in.nIter = arguments.nIter;
  in.nx = (int)floor(in.xmax/in.dx);
 
  // Solve STLS equation
  double start = omp_get_wtime();
  solveSTLS(in);
  double end = omp_get_wtime();
  printf("Solution of STLS equation complete. Elapsed time: %f seconds\n", end - start);


  return 0;

}

