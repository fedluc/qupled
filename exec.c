#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include "read_input.h"
#include "solvers.h"
#include "restart.h"
#include "dynamic_stls.h"
#include "dynamic_qstls.h"
#include "dynamic_qstls_iet.h"

void run_static_mode(input in);
void run_dynamic_mode(input in);
void run_guess_mode(input in);

int main (int argc, char **argv){

  input in;
  
  // Get input data
  get_input(argc, argv, &in);

  // Start timing
  double tic = omp_get_wtime();

  // Set number of threads for parallel calculations
  omp_set_num_threads(in.nThreads);

  // Select working mode for the code
  if (strcmp(in.mode, "static") == 0) {
    run_static_mode(in);
  }
  else if (strcmp(in.mode, "dynamic") == 0) {
    run_dynamic_mode(in);
  }
  else if (strcmp(in.mode, "guess") == 0) {
    run_guess_mode(in);
  }
  else {

    fprintf(stderr, "Error: Unknown working mode %s"
	    " Choose between: static, dynamic and guess\n"
	    , in.mode);
    exit(EXIT_FAILURE);
    
  }

  // End timing
  double toc = omp_get_wtime();

  // Print conclusion message on screen
  printf("Calculations complete. Elapsed time: %f seconds\n",
	 toc - tic);

  return 0;

}

// Run the code in "static" mode
void run_static_mode(input in){

  if (strcmp(in.theory, "STLS") == 0) {
    solve_stls(in, true);
  }
  else if (strcmp(in.theory, "VS-STLS") == 0) {
    solve_vs_stls(in, true);
  }
  else if (strcmp(in.theory, "STLS-IET-HNC") == 0 ||
	     strcmp(in.theory, "STLS-IET-IOI") == 0 ||
	   strcmp(in.theory, "STLS-IET-LCT") == 0) {
    solve_stls_iet(in,true);
  }
  else if (strcmp(in.theory, "QSTLS") == 0) {
    solve_qstls(in, true);
  }
  else if (strcmp(in.theory, "QSTLS-IET-HNC") == 0 ||
	   strcmp(in.theory, "QSTLS-IET-IOI") == 0 ||
	   strcmp(in.theory, "QSTLS-IET-LCT") == 0) {
    solve_qstls_iet(in,true);
  }
  else {
    fprintf(stderr, "Error: %s is an unknown theory to be solved. "
	    "Choose between: STLS, VS-STLS, STLS-IET-HNC,"
	    " STLS-IET-IOI, STLS-IET-LCT, QSTLS, QSTLS-IET-HNC,"
	    "QSTLS-IET-IOI and QSTLS-IET-LCT\n", in.theory);
    exit(EXIT_FAILURE);
  }
  
}

// Run the code in "dynamic" mode
void run_dynamic_mode(input in){

  if (strcmp(in.theory, "STLS") == 0 ||
      strcmp(in.theory, "VS-STLS") == 0 ||
      strcmp(in.theory, "STLS-IET-HNC") == 0 ||
      strcmp(in.theory, "STLS-IET-IOI") == 0 ||
      strcmp(in.theory, "STLS-IET-LCT") == 0) {
    compute_dynamic_stls(in, true);
  }
  else if (strcmp(in.theory, "QSTLS") == 0){
    compute_dynamic_qstls(in, true);
  }
  else if (strcmp(in.theory, "QSTLS-IET-HNC") == 0 ||
	   strcmp(in.theory, "QSTLS-IET-IOI") == 0 ||
	   strcmp(in.theory, "QSTLS-IET-LCT") == 0) {
    compute_dynamic_qstls_iet(in,true);
  }
  else {
    fprintf(stderr, "Error: %s is an unknown theory to be solved. "
	    "Choose between: STLS, VS-STLS, STLS-IET-HNC,"
	    " STLS-IET-IOI, STLS-IET-LCT, QSTLS, QSTLS-IET-HNC,"
	    "QSTLS-IET-IOI and QSTLS-IET-LCT\n", in.theory);
    exit(EXIT_FAILURE);
  }
  
}


// Run the code in guess mode
void run_guess_mode(input in){

  create_restart(in);
      
}
