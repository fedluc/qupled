#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "utils.h"

// -------------------------------------------------------------------
// LOCAL FUNCTIONS
// -------------------------------------------------------------------

// Internal energy
static double uex(double yy, void* pp);

// Radial distribution function
static double  xssf(double xx, void *pp);

// -------------------------------------------------------------------
// LOCAL DATA STRUCTURES
// -------------------------------------------------------------------

// Parameters for the integration in the internal energy expression
struct uex_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};

// Parameters for the integration in the radial distribution function
struct rdf_params {

  double cutoff;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  
};

// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A TWO-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx2(int xx, int yy, int x_size) {
  return (yy * x_size) + xx;
}


// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A THREE-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx3(int xx, int yy, int zz,
         int x_size, int y_size) {
  return (zz * x_size * y_size) + (yy * x_size) + xx;
}


// -------------------------------------------------------------------
// FUNCTION USED TO GET THE SIGN OF A NUMBER
// -------------------------------------------------------------------

int get_sign(double num) {

  if (num < 0)
    return -1;
  else
    return 1;

}  


// -------------------------------------------------------------------
// FUNCTIONS USED TO READ DATA FROM TEXT FILES
// -------------------------------------------------------------------

void get_data_format_from_text(char * file_name, int *n_lines,
			       int *n_columns){

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
    fprintf(stderr,"Error while opening text file %s\n", file_name);
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
      fprintf(stderr,"Error while reading line %d of file %s." 
	      " Only %ld elements where read\n",
	      *n_lines, file_name, num_el);
      exit(EXIT_FAILURE);
    }
    n_columns_check = 0;
    
  }

  // Close file 
  fclose(fid);
  
}

void get_data_from_text(char *file_name, int n_lines, int n_columns,
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
    fprintf(stderr,"Error while opening text file %s\n", file_name);
    exit(EXIT_FAILURE);
  }

  // Read file line by line
  for (int ii=0; ii<n_lines; ii++) {

    // Get line
    num_el = getline(&line, &len, fid);
    if (num_el == -1){
      fprintf(stderr,"Error while reading text file %s\n", file_name);
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

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

double compute_internal_energy(double *SS, double *xx,  input in) {

  double err;
  size_t neval;
  double ie;
  double lambda = pow(4.0/(9.0*M_PI), 1.0/3.0);  
  
  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &uex;

  // Internal energy
  struct uex_params uexp = {ssf_sp_ptr, ssf_acc_ptr};
  ff_int.params = &uexp;  
  gsl_integration_cquad(&ff_int,
			xx[0], xx[in.nx-1],
			0.0, QUAD_REL_ERR,
			wsp,
			&ie, &err, &neval);
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(ssf_sp_ptr);
  gsl_interp_accel_free(ssf_acc_ptr);

  // Output
  return ie/(M_PI*in.rs*lambda);

}

double uex(double yy, void* pp) {

  struct uex_params* params = (struct uex_params*)pp;
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);

  return gsl_spline_eval(ssf_sp_ptr, yy, ssf_acc_ptr) - 1;
}

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE RADIAL DISTRIBUTION FUNCTION
// -------------------------------------------------------------------

void compute_rdf(double *gg, double *rr, double *SS, double *xx, input in){

  double err;

  // Declare accelerator and spline objects
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  ssf_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nx);
  ssf_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(ssf_sp_ptr, xx, SS, in.nx);
  
  // Integration workspace
  gsl_integration_workspace *wsp 
    = gsl_integration_workspace_alloc(1000);
  gsl_integration_workspace *wspc 
    = gsl_integration_workspace_alloc(1000);
  gsl_integration_qawo_table *qtab
    = gsl_integration_qawo_table_alloc(0.0, 1.0, GSL_INTEG_SINE, 1000);
  
  // Integration function
  gsl_function ff_int;
  struct rdf_params rdfp = {xx[in.nx-1], ssf_sp_ptr, ssf_acc_ptr};
  ff_int.function = &xssf;
  ff_int.params = &rdfp;

  // Real space grid
  for (int ii = 0; ii < in.nx; ii++){
    rr[ii] = 0.01 + in.dx*ii;
  }
  
  // Radial distribution function
  for (int ii = 0; ii < in.nx; ii++) {

    // Set wave-vector (divide xx[ii] by ll to convert to Wigner-Seitz units)
    gsl_integration_qawo_table_set(qtab, rr[ii], 1.0, GSL_INTEG_SINE);

    // Fourier transform
    gsl_integration_qawf(&ff_int,
    			 0.0, QUAD_REL_ERR,
			 1000,
    			 wsp, wspc,
    			 qtab,
    			 &gg[ii], &err);
    gg[ii] = 1 + 3.0/(2.0*rr[ii])*gg[ii];

  }

  // Free memory
  gsl_integration_workspace_free(wsp);
  gsl_integration_workspace_free(wspc); 
  gsl_integration_qawo_table_free(qtab);
  gsl_spline_free(ssf_sp_ptr);
  gsl_interp_accel_free(ssf_acc_ptr);
 

}

double  xssf(double xx, void *pp){
  
  struct rdf_params* params = (struct rdf_params*)pp;
  double cutoff = (params->cutoff);
  gsl_spline* ssf_sp_ptr = (params->ssf_sp_ptr);
  gsl_interp_accel* ssf_acc_ptr = (params->ssf_acc_ptr);

  if (xx < cutoff)
    return xx*(gsl_spline_eval(ssf_sp_ptr, xx, ssf_acc_ptr) - 1);
  else
    return 0.0;
}
