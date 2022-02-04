#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "restart.h"
#include "stls.h"
#include "qstls.h"

// -------------------------------------------------------------------
// FUNCTION USED TO WRITE BINARY FILES FOR GUESS (OR RESTART) STARTING
// FROM TEXT FILES OBTAINED FROM A SIMULATION
// -------------------------------------------------------------------

void create_restart(input in){

  // Variables
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
  
  // Get format of the data stored in the text files
  get_restart_data_format(in.guess_file1, &n_lines_file1, &n_columns_file1);
  get_restart_data_format(in.guess_file2, &n_lines_file2, &n_columns_file2);

  // Update input structure based on the content of the text files
  set_nx_nl(n_lines_file1, n_lines_file2, n_columns_file1, n_columns_file2, &in);

  // Allocate stls arrays
  alloc_stls_arrays(in, &xx, &phi, &GG, &GG_new, &SS, &SSHF);

  // Allocate additional qstsl arrays if necessary
  if(in.theory_id > 5) {
    psi = malloc( sizeof(double) * in.nx * in.nl);
  }

  // Get restart data
  if(in.theory_id <= 5) {

    // Static structure factor
    get_restart_data(in.guess_file1, n_lines_file1, n_columns_file1, SS, xx, &in);
    // Static local field correction
    get_restart_data(in.guess_file2, n_lines_file2, n_columns_file2, GG, xx, &in);

  }
  else{

    // Static structure factor
    get_restart_data(in.guess_file1, n_lines_file1, n_columns_file1, SS, xx, &in);
    // Auxiliary density response
    get_restart_data(in.guess_file2, n_lines_file2, n_columns_file2, psi, NULL, &in);

  }
  
  // Write restart files
  sprintf(out_name, "restart_rs%.3f_theta%.3f_%s.bin", in.rs, in.Theta, in.theory);
  if(in.theory_id <= 5) {
    write_guess_stls(SS, GG, in);
  }
  else{
    write_guess_qstls(SS, psi, in);
  }

  
  // Write concluding message on screen
  printf("Binary files successfully created. See infomation below\n");
  printf("Theory: %s\n", in.theory);
  printf("Quantum coupling parameter: %f\n", in.rs);
  printf("Quantum degeneracy parameter: %f\n", in.Theta);
  printf("Number of grid points: %d\n", in.nx);
  printf("Resolution for wave-vector grid: %f\n", in.dx);
  printf("Cutoff for wave-vector grid: %f\n", in.xmax);
  if(in.theory_id <= 5) {
    printf("Source file for the static structure factor: %s\n", in.guess_file1);
    printf("Source file for the static local field correction: %s\n", in.guess_file2);
  }
  else{
    printf("Number of Matsubara frequencies: %d\n", in.nl);
    printf("Source file for the static structure factor: %s\n", in.guess_file1);
    printf("Source file for the static local field correction: %s\n", in.guess_file2);
  }
  printf("Output file: %s\n", out_name); 
  
  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
  free(psi);
  
}
  

// -------------------------------------------------------------------
// FUNCTION USED TO INFER THE FORMAT OF THE TEXT FILES
// -------------------------------------------------------------------
void get_restart_data_format(char * file_name, int *n_lines, int *n_columns){

  // Variables
  FILE *fid;
  char * line = NULL;
  char * value = NULL;
  size_t len = 0;
  ssize_t num_el;
  int n_columns_check;
  
  // Open file
  fid = NULL;
  fid = fopen(file_name, "r");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file %s to construct the restart\n", file_name);
    exit(EXIT_FAILURE);
  }

  // Get number of lines
  *n_lines = 0;
  *n_columns = 0;
  n_columns_check = 0;
  while ((num_el = getline(&line, &len, fid)) != -1) {

    // Update line counter
    *n_lines += 1;

    // Get number of columns and check that format remains consistent
    value = strtok(line, " \n");
    while(value != NULL){
      if(*n_lines == 1) *n_columns += 1;
      else n_columns_check += 1;
      value = strtok(NULL, " \n");
    }
    if (*n_lines > 1 && n_columns_check != *n_columns){
      fprintf(stderr,"Error while reading line %d of file %s. Only %ld elements where read\n",
	      *n_lines, file_name, num_el);
      exit(EXIT_FAILURE);
    }
    n_columns_check = 0;
    
  }

  // Close file 
  fclose(fid);
  
}


// -------------------------------------------------------------------
// FUNCTION USED TO SET THE PARAMETERS FOR THE GRID SIZE AND FOR THE
// NUMBER OF MATSUBARA FREQUENCIES
// -------------------------------------------------------------------
void set_nx_nl(int nl1, int nl2, int nc1, int nc2, input *in){

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

  if (nc2>2 && in->theory_id <= 5){
    fprintf(stderr, "Unexpected data format for %s theory files. "
	   "Only two columns are expected, but %d where detected\n", in->theory, nc2);
    exit(EXIT_FAILURE);
  }
  
}


// -------------------------------------------------------------------
// FUNCTION USED TO READ THE TEXT FILES
// -------------------------------------------------------------------
void get_restart_data(char *file_name, int n_lines, int n_columns,
		      double *data, double *xx, input *in){

  // Variables
  FILE *fid;
  char * line = NULL;
  char * value = NULL;
  size_t len = 0;
  ssize_t num_el;

  // Open file
  fid = NULL;
  fid = fopen(file_name, "r");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening file %s to construct the restart\n", file_name);
    exit(EXIT_FAILURE);
  }

  // Read file line by line
  for (int ii=0; ii<n_lines; ii++) {

    // Get line
    num_el = getline(&line, &len, fid);
    if (num_el == -1){
      fprintf(stderr,"Error while reading file %s during data extraction\n", file_name);
      exit(EXIT_FAILURE);
    }
    
    // Extract information from the first column
    value = strtok(line, " \n");
    if (xx == NULL) {
      // First column is data
      data[idx2(ii,0,n_lines)] = atof(value);
    }
    else{
      // First column is the grid
      xx[ii] = atof(value);
    }

    // Extract information from the remaining columns
    for (int jj=1; jj<n_columns; jj++) {
      value = strtok(NULL, " \n");
      if (xx == NULL) {
	data[idx2(ii,jj,n_lines)] = atof(value);
      }
      else{
	data[idx2(ii,jj-1,n_lines)] = atof(value);
      }
    }
    
  }

  // Close file 
  fclose(fid);

  // Update input structure
  if (xx != NULL) {
    in->xmax = xx[in->nx - 1];
    in->dx = xx[1] - xx[0];
  }
    
}



