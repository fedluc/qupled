#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "restart.h"
#include "stls.h"
#include "qstls.h"

// -------------------------------------------------------------------
// LOCAL FUNCTIONS
// -------------------------------------------------------------------

// Set wave-vector grid size and Matsubara frequencies
static void set_nx_nl(int nl1, int nl2, int nc1, int nc2,
		      bool is_qstls, input *in);

// -------------------------------------------------------------------
// FUNCTION USED TO WRITE BINARY FILES FOR GUESS (OR RESTART) STARTING
// FROM TEXT FILES OBTAINED FROM A SIMULATION
// -------------------------------------------------------------------

void create_restart(input in){

  // Variables
  bool is_qstls;
  int n_lines_file1;
  int n_lines_file2;
  int n_columns_file1;
  int n_columns_file2;
  double *xx = NULL; 
  double *phi = NULL;
  double *GG = NULL;
  double *GG_new = NULL;
  double *SS = NULL;
  double *SSHF = NULL;
  double *psi = NULL;
  char out_name[100];

  // Check for which theory the restart file must be prepared
  if (strstr(in.theory, "QSTLS") != NULL) is_qstls = true;
  else is_qstls = false;
  
  // Get format of the data stored in the text files
  get_data_format_from_text(in.guess_file1, &n_lines_file1, &n_columns_file1);
  get_data_format_from_text(in.guess_file2, &n_lines_file2, &n_columns_file2);

  // Update input structure based on the content of the text files
  set_nx_nl(n_lines_file1, n_lines_file2, n_columns_file1,
	    n_columns_file2, is_qstls, &in);

  // Allocate stls arrays
  alloc_stls_arrays(in, &xx, &phi, &GG, &GG_new, &SS, &SSHF);

  // Allocate additional qstsl arrays if necessary
  if(is_qstls) {
    psi = malloc( sizeof(double) * in.nx * in.nl);
    if (psi == NULL) {
      fprintf(stderr, "Failed to allocate memory for the auxiliary density response\n");
      exit(EXIT_FAILURE);
    }
  }

  // Get restart data
  if(is_qstls) {

    // Static structure factor
    get_data_from_text(in.guess_file1, n_lines_file1, n_columns_file1, SS, xx, &in);
    // Auxiliary density response
    get_data_from_text(in.guess_file2, n_lines_file2, n_columns_file2, psi, NULL, &in);


  }
  else{
    
    // Static structure factor
    get_data_from_text(in.guess_file1, n_lines_file1, n_columns_file1, SS, xx, &in);
    // Static local field correction
    get_data_from_text(in.guess_file2, n_lines_file2, n_columns_file2, GG, xx, &in);

  }
  
  // Write restart files
  sprintf(out_name, "restart_rs%.3f_theta%.3f_%s.bin", in.rs, in.Theta, in.theory);
  if(is_qstls) {
    write_guess_qstls(SS, psi, in);
  }
  else{
    write_guess_stls(SS, GG, in);
  }

  
  // Write concluding message on screen
  printf("Binary files successfully created. See infomation below\n");
  printf("Theory: %s\n", in.theory);
  printf("Quantum coupling parameter: %f\n", in.rs);
  printf("Quantum degeneracy parameter: %f\n", in.Theta);
  printf("Number of grid points: %d\n", in.nx);
  printf("Resolution for wave-vector grid: %f\n", in.dx);
  printf("Cutoff for wave-vector grid: %f\n", in.xmax);
  if(is_qstls) {
    printf("Number of Matsubara frequencies: %d\n", in.nl);
    printf("Source file for the static structure factor: %s\n", in.guess_file1);
    printf("Source file for the static local field correction: %s\n", in.guess_file2);
  }
  else{
    printf("Source file for the static structure factor: %s\n", in.guess_file1);
    printf("Source file for the static local field correction: %s\n", in.guess_file2); 
  }
  printf("Output file: %s\n", out_name); 
  
  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
  if (is_qstls) free(psi);
  
}
  

// -------------------------------------------------------------------
// FUNCTION USED TO SET THE PARAMETERS FOR THE GRID SIZE AND FOR THE
// NUMBER OF MATSUBARA FREQUENCIES
// -------------------------------------------------------------------
void set_nx_nl(int nl1, int nl2, int nc1, int nc2, bool is_qstls, input *in){

  if (nl1 == nl2) {
    in->nx = nl2;
    in->nl = nc2; 
  }
  else {
    fprintf(stderr,"The files used to construct the restart are inconsistent\n");
    fprintf(stderr, "File: %s\n %d lines, %d columns\n", in->stls_guess_file, nl1, nc1);
    fprintf(stderr, "File: %s\n %d lines, %d columns\n", in->stls_guess_file, nl2, nc2);
    exit(EXIT_FAILURE);
  }

  if (nc2>2 && !is_qstls){
    fprintf(stderr, "Unexpected data format for %s theory files. "
	   "Only two columns are expected, but %d where detected\n", in->theory, nc2);
    exit(EXIT_FAILURE);
  }
  
}


