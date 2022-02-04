#include <stdio.h>
#include <omp.h>
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
  }
  // Compute new solution
  else{

    // Set number of threads for parallel calculations
    omp_set_num_threads(in.nThreads);

    // Solve theory
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
    
  }
    

  // End timing
  double toc = omp_get_wtime();

  // Print conclusion message on screen
  printf("Solution complete. Elapsed time: %f seconds\n", toc - tic);

  return 0;

}
