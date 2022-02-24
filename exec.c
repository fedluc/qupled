#include <stdio.h>
#include <omp.h>
#include <string.h>
#include "read_input.h"
#include "solvers.h"
#include "restart.h"


int main (int argc, char **argv){

  input in;
  
  // Get input data
  get_input(argc, argv, &in);

  // Start timing
  double tic = omp_get_wtime();

  // Write guess (or restart) files
  if (in.guess_write) {
    create_restart(in);
    return 0;
  }
  
  // Set number of threads for parallel calculations
  omp_set_num_threads(in.nThreads);

  // Solve theory specified in input
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
    printf("Error: %s is an unknown theory to be solved. "
  	   "Choose between: STLS, VS-STLS, STLS-IET-HNC,"
  	   " STLS-IET-IOI, STLS-IET-LCT, QSTLS, QSTLS-IET-HNC,"
  	   "QSTLS-IET-IOI and QSTLS-IET-LCT\n", in.theory);
    return 1;
  }
  
  // End timing
  double toc = omp_get_wtime();

  // Print conclusion message on screen
  printf("Solution complete. Elapsed time: %f seconds\n", toc - tic);

  return 0;

}
