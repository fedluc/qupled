#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "utils.h"

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
// FUNCTIONS USED TO COMPUTE THE INTERNAL ENERGY
// -------------------------------------------------------------------

struct uex_params {

  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;

};

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
// FUNCTIONS USED TO COMPUTE THE FREE ENERGY
// -------------------------------------------------------------------

struct fex_params {

  gsl_spline *rsu_sp_ptr;
  gsl_interp_accel *rsu_acc_ptr;

};

double compute_free_energy(double *rsu, double *rsp, input in) {

  double err;
  size_t neval;
  double fre;

  // Declare accelerator and spline objects
  gsl_spline *rsu_sp_ptr;
  gsl_interp_accel *rsu_acc_ptr;
  
  // Allocate the accelerator and the spline objects
  rsu_sp_ptr = gsl_spline_alloc(gsl_interp_cspline, in.nrs);
  rsu_acc_ptr = gsl_interp_accel_alloc();
  
  // Initialize the spline
  gsl_spline_init(rsu_sp_ptr, rsp, rsu, in.nrs);

  // Integration workspace
  gsl_integration_cquad_workspace *wsp
    = gsl_integration_cquad_workspace_alloc(100);

  // Integration function
  gsl_function ff_int;
  ff_int.function = &fex;

  // Internal energy
  struct fex_params fexp = {rsu_sp_ptr, rsu_acc_ptr};
  ff_int.params = &fexp;
  gsl_integration_cquad(&ff_int,
  			rsp[0], in.rs,
  			0.0, QUAD_REL_ERR,
  			wsp,
  			&fre, &err, &neval);
  
  // Free memory
  gsl_integration_cquad_workspace_free(wsp);
  gsl_spline_free(rsu_sp_ptr);
  gsl_interp_accel_free(rsu_acc_ptr);

  // Output
  return fre/(in.rs*in.rs);

}

double fex(double rs, void* pp) {

  struct fex_params* params = (struct fex_params*)pp;
  gsl_spline* rsu_sp_ptr = (params->rsu_sp_ptr);
  gsl_interp_accel* rsu_acc_ptr = (params->rsu_acc_ptr);

  return gsl_spline_eval(rsu_sp_ptr, rs, rsu_acc_ptr);
  
}



// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE RADIAL DISTRIBUTION FUNCTION
// -------------------------------------------------------------------

struct rdf_params {

  double cutoff;
  gsl_spline *ssf_sp_ptr;
  gsl_interp_accel *ssf_acc_ptr;
  
};

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
