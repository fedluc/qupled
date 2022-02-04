#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <argp.h>
#include "read_input.h"

// ----------------------------------------
// CONSTANTS AND DATA STRUCTURES
// ----------------------------------------

static bool debug_input;

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
#define ARGUMENT_DEBUG_SHORT 0x95
#define ARGUMENT_GUESS_WRITE_SHORT 0x96
#define ARGUMENT_GUESS_FILES_SHORT 0x97


// Optional arguments
static struct argp_option options[] = {
				       
  {"Theta", ARGUMENT_THETA_SHORT, "1.0", 0,
   "Quantum degeneracy parameter"},
  
  {"rs", ARGUMENT_RS_SHORT, "1.0", 0,
   "Quantum coupling parameter"},
  
  {"xmax", ARGUMENT_XMAX_SHORT, "20.0", 0,
   "Cutoff for wave-vector grid"},
  
  {"dx", ARGUMENT_DX_SHORT, "0.1", 0,
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

  {"theory", ARGUMENT_THEORY_SHORT, "STLS",0,
   "Scheme to be solved"},

  {"omp", ARGUMENT_OMP_SHORT, "1",0,
   "Number of omp threads to use in the solution"},

  {"debug-input", ARGUMENT_DEBUG_SHORT, "0",0,
   "Print content of the input structure on screen (0 = off, 1 = on)"},

  {"guess-write", ARGUMENT_GUESS_WRITE_SHORT, "0",0,
   "Write binary restart files from text files (0 = off, 1 = on)"},

  {"guess-files", ARGUMENT_GUESS_FILES_SHORT, "NO_FILE,NO_FILE",0,
   "Text files to write binary restart files"},

  { 0 }
  
};


// Command line parser
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{

  char *value;
  input *in = state->input;
  
  switch (key)
    {

    case ARGUMENT_DX_SHORT:
      in->dx = atof(arg);
      break;
      
    case  ARGUMENT_MIN_ERR_SHORT:
      in->err_min_iter = atof(arg);
      break;
      
    case  ARGUMENT_STLS_GUESS_SHORT:
      in->stls_guess_file = arg;
      break;
      
    case  ARGUMENT_QSTLS_GUESS_SHORT:
      in->qstls_guess_file = arg;
      break;
      
    case  ARGUMENT_QSTLS_FIXED_SHORT:
      in->qstls_fixed_file = arg;
      break;
      
    case  ARGUMENT_QSTLS_IET_FIXED_SHORT:
      in->qstls_iet_fixed_file = arg;
      break;
      
    case  ARGUMENT_MU_GUESS_SHORT:
      value = strtok(arg, ",");
      if(value != NULL ) {
	in->mu_lo = atof(value);
      }
      else exit(EXIT_FAILURE);
      value = strtok(NULL, ",");
      if(value != NULL ) {
	in->mu_hi = atof(value);
      }
      else exit(EXIT_FAILURE);
      break;
      
    case ARGUMENT_ITER_SHORT:
      in->nIter = atoi(arg);
      break;
      
    case ARGUMENT_NL_SHORT:
      in->nl = atoi(arg);
      break;
      
    case ARGUMENT_MIX_SHORT:
      in->a_mix = atof(arg);
      break;
      
    case ARGUMENT_OMP_SHORT:
      in->nThreads = atoi(arg);
      break;
      
    case ARGUMENT_RS_SHORT:
      in->rs = atof(arg);
      break;
      
    case ARGUMENT_THETA_SHORT:
      in->Theta = atof(arg);
      break;
      
    case ARGUMENT_THEORY_SHORT:
      in->theory = arg;
      break;
      
    case  ARGUMENT_XMAX_SHORT:
      in->xmax = atof(arg);
      break;

    case  ARGUMENT_DEBUG_SHORT:
      debug_input = atoi(arg);
      break;

    case  ARGUMENT_GUESS_WRITE_SHORT:
      in->guess_write = atoi(arg);
      break;  

    case  ARGUMENT_GUESS_FILES_SHORT:
      value = strtok(arg, ",");
      if(value != NULL ) {
	in->guess_file1 = value;
      }
      else exit(EXIT_FAILURE);
      value = strtok(NULL, ",");
      if(value != NULL ) {
	in->guess_file2 = value;
      }
      else exit(EXIT_FAILURE);
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
// FUNCTION TO READ INPUT DATA
// ----------------------------------------
void get_input(int argc, char **argv, input *in){

  // Default values for optional arguments
  set_default_parse_opt(in);

  // Parse command line
  argp_parse(&argp, argc, argv, 0, 0, in);

  // Get number of grid points
  get_nx(in);
  
  // Assign a numerical id to the theory given in input
  get_theory_id(in);

  // Debug input
  if(debug_input) print_input(in);
  
}


// -------------------------------------------------
// FUNCTION TO ASSIGN DEFAULT VALUES TO PARSER DATA
// -------------------------------------------------
void set_default_parse_opt(input *in){

  debug_input = false;
  in->stls_guess_file = "NO_FILE"; // File with initial guess for STLS and STLS-IET schemes
  in->qstls_guess_file = "NO_FILE"; // File with initial guess for QSTLS and QSTLS-IET schemes
  in->qstls_fixed_file = "NO_FILE"; // File with fixed component of the density response for the QSTLS scheme
  in->qstls_iet_fixed_file = "NO_FILE"; // File with fixed component of the density response for the QSTLS-IET scheme
  in->Theta = 1.0; // Quantum degeneracy parameter
  in->rs = 1.0; // Quantum coupling parameter
  in->dx = 0.1; // Wave-vector grid resolution
  in->err_min_iter = 1e-5; // Minimum error for convergence in the iterative procedure
  in->a_mix = 0.1; // Mixing parameter for iterative procedure
  in->mu_lo = -10; // Initial guess for chemical potential (low bound)
  in->mu_hi = 10; // Initial guess for chemical potential (high bound)
  in->xmax = 20; // Cutoff for wave-vector grid
  in->nl = 128; // Number of Matsubara frequencies
  in->nIter = 1000; // Number of iterations
  in->theory = "STLS"; // Theory to solve
  in->nThreads = 1; // Number of OMP threads to use in the solution
  in->guess_write = 0; // Flag to run in "write guess" mode
  in->guess_file1 = "NO_FILE"; // File of the first file used to construct the guess (static structure factor)
  in->guess_file2 = "NO_FILE"; // File of the second file used to construct the guess (static local field correction or auxiliary density response)
    
  
}

// ------------------------------------------------------------------
// FUNCTION TO ASSIGN A NUMERICAL ID TO THE THEORY SPECIFIED IN INPUT
// ------------------------------------------------------------------
void get_theory_id(input *in){

  if (strcmp(in->theory, "STLS") == 0) in->theory_id = 1;
  else if (strcmp(in->theory, "STLS-IET-HNC") == 0) in->theory_id = 2;
  else if (strcmp(in->theory, "STLS-IET-IOI") == 0) in->theory_id = 3;
  else if (strcmp(in->theory, "STLS-IET-LCT") == 0) in->theory_id = 4;
  else if (strcmp(in->theory, "STLS-RIET-LCT") == 0) in->theory_id = 5;
  else if (strcmp(in->theory, "QSTLS") == 0) in->theory_id = 6;
  else if (strcmp(in->theory, "QSTLS-IET-HNC") == 0) in->theory_id = 7;
  else if (strcmp(in->theory, "QSTLS-IET-IOI") == 0) in->theory_id = 8;
  else if (strcmp(in->theory, "QSTLS-IET-LCT") == 0) in->theory_id = 9;
  else if (strcmp(in->theory, "QSTLS-RIET-LCT") == 0) in->theory_id = 10;
  else {
    printf("Error: unknown theory to be solved. " 
	   "Choose between: STLS, STLS-IET-HNC, STLS-IET-IOI,"
	   "STLS-IET-LCT, STLS-RIET-LCT, QSTLS, QSTLS-IET-HNC,"
	   "QSTLS-IET-IOI, QSTLS-IET-LCT and QSTLS-RIET-LCT\n");
    exit(EXIT_FAILURE);
  }
}

// ------------------------------------------------
// FUNCTION TO COMPUTE THE NUMBER OF GRID POINTS
// ------------------------------------------------
void get_nx(input *in){
  in->nx = (int)floor(in->xmax/in->dx);
}


// ------------------------------------------------
// FUNCTION TO DEBUG THE INPUT
// ------------------------------------------------
void print_input(input *in){
  
  printf("------ Input parameters -------------\n");
  printf("File for initial guess (STLS): %s\n", in->stls_guess_file);
  printf("File for initial guess (qSTLS): %s\n", in->qstls_guess_file);
  printf("File for fixed component (qSTLS): %s\n", in->qstls_fixed_file);
  printf("File for fixed component (qSTLS-IET): %s\n", in->qstls_iet_fixed_file);
  printf("Theory: %s\n", in->theory);
  printf("Quantum degeneracy parameter: %f\n", in->Theta);
  printf("Quantum coupling parameter: %f\n", in->rs);
  printf("Wave-vector resolutions: %f\n", in->dx);
  printf("Error for convergence: %.5e\n", in->err_min_iter);
  printf("Mixing parameter: %f\n", in->a_mix);
  printf("Chemical potential (low and high bound): %f %f\n", 
	 in->mu_lo, in->mu_hi);
  printf("Wave-vector cutoff: %f\n", in->xmax);
  printf("Number of Matsubara frequencies: %d\n", in->nl);
  printf("Number of grid points: %d\n", in->nIter);  
  printf("Maximum number of iterations: %d\n", in->nIter);
  printf("Number of threads: %d\n", in->nThreads);
  printf("Theory ID: %d\n", in->theory_id);
  printf("Write-guess mode: %d\n", in->guess_write);
  printf("Guess file 1: %s\n", in->guess_file1);
  printf("Guess file 2: %s\n", in->guess_file2);
  printf("-------------------------------------\n");
  
}
