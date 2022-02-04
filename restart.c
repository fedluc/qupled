#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "restart.h"
#include "stls.h"
#include "qstls.h"

// -------------------------------------------------------------------
// FUNCTION USED TO ...
// -------------------------------------------------------------------

// 
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
  
  // Get format of the data stored in the text files
  get_restart_data_format(in.stls_guess_file, &n_lines_file1, &n_columns_file1);
  get_restart_data_format(in.stls_guess_file, &n_lines_file2, &n_columns_file2);

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
    get_restart_data(in.stls_guess_file, n_lines_file1, n_columns_file1, SS, xx, &in);
    // Static local field correction
    //get_restart_data(in.stls_guess_file, n_lines_file2, n_columns_file2, GG, xx, &in);

  }
  else{

    // Static structure factor
    get_restart_data(in.stls_guess_file, n_lines_file1, n_columns_file1, SS, xx, &in);
    // Auxiliary density response
    get_restart_data(in.stls_guess_file, n_lines_file2, n_columns_file2, GG, NULL, &in);

  }
  
  // Write restart files
  if(in.theory_id <= 5) {
    write_guess_stls(SS, GG, in);
  }
  else{
    write_guess_qstls(SS, psi, in);
  }

  // Free memory
  free_stls_arrays(xx, phi, GG, GG_new, SS, SSHF);
  free(psi);
  
}
  

void get_restart_data_format(char * file_name, int *n_lines, int *n_columns){

  // Variables
  FILE *fid;
  char * line = NULL;
  char * value = NULL;
  size_t len = 0;
  ssize_t num_el;
  ssize_t num_el_ref = -1;
  
  
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
  while ((num_el = getline(&line, &len, fid)) != -1) {

    // Check that line format remains consistent
    if (num_el_ref < 0) {
      num_el_ref = num_el;
    }
    if (num_el != num_el_ref){
      fprintf(stderr,"Error while reading file %s. Inconsistent file format. \n", file_name);
      exit(EXIT_FAILURE);
    }

    // Update line counter
    *n_lines += 1;

    // Get number of columns (only done for the first line)
    if (*n_lines == 1) {
      value = strtok(line, " \n");
      while(value != NULL){
      	*n_columns += 1;
      	value = strtok(NULL, " \n");
      }
    }
    
    }

  // Close file 
  fclose(fid);
  
}



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
  
}



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



