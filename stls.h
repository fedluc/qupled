#ifndef STLS_H
#define STLS_H

#include <stdbool.h>

typedef struct {

  char *phi_file;
  char *ssf_file;
  double Theta;
  double rs;
  double dx;
  double err_min_iter;
  double a_mix;
  double mu_lo;
  double mu_hi;
  double mu;
  double xmax;
  int nl;
  int nx;
  int nIter;


} input;


void solve_stls(input in, bool verbose,
                double **xx_out, double **SS_out,
                double **SSHF_out, double **GG_out,
                double **GG_new_out, double **phi_out);

void alloc_stls_arrays(input in, double **xx, double **phi,
		       double **GG, double **GG_new,
		       double **SS, double **SSHF);

void free_stls_arrays(double *xx, bool free_xx,
                      double *phi, bool free_phi,
                      double *GG, bool free_GG,
                      double *GG_new, bool free_GG_new,
                      double *SS, bool free_SS,
                      double *SSHF, bool free_SSHF);

double compute_mu(input in);

double normalization_condition(double mu, void *pp);

void wave_vector_grid(double *xx, input in);

void compute_phi(double *phi, double *xx, input in, bool verbose);

void compute_phil(double *phil, double *xx, int ll, input in);

double phixl(double yy, double xx, int ll, input in);

double phix0(double yy, double xx, input in);

int idx2(int xx, int yy, int x_size);

void compute_ssfHF(double *SS,  double *xx, input in);

double ssfHF(double yy, double xx, input in);


void compute_ssf(double *SS, double *SSHF,
                 double *GG, double *phi, 
		 double *xx, input in);

void compute_slfc(double *GG, double *SS, 
		  double *xx, input in);

double slfc(double yy, double xx, double SS);

double compute_uex(double *SS, input in);

void write_text(double *SS, double *GG, 
		double *xx, input in );

void write_bin(double *phi, double *SSHF, input in);

void read_text(double *SS, double *GG, 
	       double *xx, input in);

void read_bin(input *in, double **xx, double **phi,
              double **GG, double **GG_new,
              double **SS, double **SSHF);

#endif
