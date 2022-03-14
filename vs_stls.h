#ifndef VS_STLS_H
#define VS_STLS_H

#include "read_input.h"


// -------------------------------------------------------------------
// DATA STRUCTURES TO HANDLE MORE STATE POINTS SIMULTANEOUSLY
// -------------------------------------------------------------------

#define VSS_NUMEL 9
#define VSS_IDXIN 4
#define VSS_STENCIL 3
#define VST_NUMEL 3
#define VST_IDXIN 1

typedef union {

  struct {

    double *rsmtm;
    double *rstm;
    double *rsptm;

    double *rsmt;
    double *rst;
    double *rspt;

    double *rsmtp;
    double *rstp;
    double *rsptp;
    
  };
  double *el[VSS_NUMEL];
  
} vs_struct;


typedef union {

  struct {

    double *rstm;
    double *rst;
    double *rstp;
    
  };
  double *el[VSS_NUMEL];
  
} vs_thermo;


// -------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE SIZE OF THE GRID FOR THERMODYNAMIC
// INTEGRATION
// -------------------------------------------------------------------

void get_grid_thermo_size(input *in);
  
// -------------------------------------------------------------------
// FUNCTIONS USED TO ALLOCATE AND FREE ARRAYS
// -------------------------------------------------------------------
  
void alloc_vs_stls_arrays(input in, double **rsa, double **rsu);
  
void free_vs_stls_arrays(double *rsu, double *rsa);

// -------------------------------------------------------------------
// FUNCTION USED TO INITIALIZE ARRAYS
// -------------------------------------------------------------------

void init_fixed_vs_stls_arrays(input *in, input *vs_in,
			       vs_struct xx, vs_thermo rsa,
			       vs_thermo rsu, double *fxc,
			       double *rs_co, bool verbose);

void init_tmp_vs_stls_arrays(input *vs_in, vs_struct xx,
			     vs_struct phi, vs_struct SSHF,
			     bool verbose);

// ------------------------------------------------------------------
// FUNCTION USED TO DEFINE THE WAVE-VECTOR GRID
// ------------------------------------------------------------------

void rs_grid(double *rsp, input *in);

// ---------------------------------------------------------------------
// FUNCTIONS USED TO DEFINE THE INITIAL GUESS
// ---------------------------------------------------------------------

void initial_guess_vs_stls(vs_struct xx, vs_struct SS,
			   vs_struct SSHF, vs_struct GG,
			   vs_struct GG_new, vs_struct phi,
			   input *vs_in);

// ---------------------------------------------------------------------
// FUNCTIONS USED TO PERFORM THE ITERATIONS FOR THE VS-STLS SCHEME
// ---------------------------------------------------------------------

void vs_stls_thermo_iterations(vs_struct xx, vs_thermo rsu,
			       vs_thermo rsa, double *fxc,
			       double *rs_co, input *vs_in,
			       bool verbose);

double vs_stls_thermo_err(double alpha, input *vs_in);

void vs_stls_thermo_update(double alpha, input *vs_in);

void vs_stls_struct_iterations(vs_struct SS, vs_struct SSHF,
			       vs_struct GG, vs_struct GG_new,
			       vs_struct phi, vs_struct xx,
			       input *vs_in, bool verbose);

double vs_stls_struct_err(vs_struct GG, vs_struct GG_new,
			  input *vs_in);

void vs_stls_struct_update(vs_struct GG, vs_struct GG_new,
			   input *vs_in);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE CHEMICAL POTENTIAL
// -------------------------------------------------------------------

void compute_chemical_potential_vs_stls(input *vs_in);

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE STATIC STRUCTURE FACTOR
// -------------------------------------------------------------------

void compute_ssf_vs_stls(vs_struct SS, vs_struct SSHF, vs_struct GG,
			 vs_struct phi, vs_struct xx, input *vs_in);

void compute_ssf_HF_vs_stls(vs_struct SS, vs_struct xx, input *vs_in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE STATIC LOCAL FIELD CORRECTION
// -------------------------------------------------------------------

void compute_slfc_vs_stls(vs_struct GG, vs_struct SS,
			  vs_struct xx, input *vs_in);

void slfc_vs_stls_stls(vs_struct GG, vs_struct SS,
		       vs_struct GG_stls, vs_struct xx,
		       input *vs_in);

void slfc_vs_stls_dx(vs_struct GG, vs_struct GG_stls,
		     vs_struct xx, input *vs_in);

void slfc_vs_stls_drs(vs_struct GG, vs_struct GG_stls,
		      input *vs_in);

void slfc_vs_stls_dt(vs_struct GG, vs_struct GG_stls,
		     input *vs_in);

double slfc_vs_stls_drs_centered(vs_struct GG, int ii, int jj,
				 input *vs_in);

double slfc_vs_stls_drs_forward(vs_struct GG, int ii, int jj,
				input *vs_in);

double slfc_vs_stls_drs_backward(vs_struct GG, int ii, int jj,
				 input *vs_in);

double slfc_vs_stls_dt_centered(vs_struct GG, int ii, int jj,
				input *vs_in);

double slfc_vs_stls_dt_forward(vs_struct GG, int ii, int jj,
			       input *vs_in);

double slfc_vs_stls_dt_backward(vs_struct GG, int ii, int jj,
				input *vs_in);

// -----------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE NORMALIZED IDEAL LINDHARD DENSITY RESPONSE
// -----------------------------------------------------------------------

void compute_idr_vs_stls(vs_struct phi, vs_struct xx, input *vs_in);
 
// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE PARAMETER FOR THE CSR RULE
// -------------------------------------------------------------------

double compute_alpha(vs_struct xx, vs_thermo rsu,
		     vs_thermo rsa, double *fxc,
		     double *rs_co, input *vs_in);


// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE INTEGRAND FOR THE FREE ENERGY
// -------------------------------------------------------------------

void compute_rsu(vs_struct xx, vs_thermo rsu,
		 vs_thermo rsa, double *rs_co,
		 input *vs_in,
		 bool verbose);

void compute_rsu_blocks(vs_struct SS, vs_struct SSHF,
			vs_struct GG, vs_struct GG_new,
			vs_struct phi, vs_struct xx,
			vs_thermo rsu, vs_thermo rsa,
			double *rs_co, int start,
			int end, int step, input *vs_in,
			bool compute_guess, bool verbose);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FREE ENERGY
// -------------------------------------------------------------------

double compute_free_energy(double *rsa, double *rsu,
			   double rs_min, double rs_max,
			   input in);

double fxc(double rs, void* pp);

void compute_free_energy_fixed(vs_thermo rsa, vs_thermo rsu,
			       double *fxc, double *rs_co,
			       input *vs_in);

// -------------------------------------------------------------------
// FUNCTIONS FOR OUTPUT AND INPUT
// -------------------------------------------------------------------

void write_text_vs_stls(double *SS, double *GG, double *phi,
			double *SSHF, double *xx, double *rsu,
			double *rsa, input in);

void write_text_fxc(double *rsu, double *rsp, input in);

void write_text_alpha_CSR(input in);

void read_guess_vs_stls(vs_struct SS, vs_struct GG, input *vs_in);

void write_thermo_vs_stls(vs_thermo rsa, vs_thermo rsu, input *vs_in);

void read_thermo_vs_stls(vs_thermo *rsa, vs_thermo *rsu,
			 int *nrs, input *vs_in);

void check_thermo_vs_stls(double dt, double Theta, input in,
			  size_t it_read, size_t it_expected, FILE *fid,
			  bool check_grid, bool check_items, bool check_eof);
#endif

