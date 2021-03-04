#ifndef STLS_H
#define STLS_H

typedef struct {

  char *phi_file;
  char *ssf_file;
  double Theta;
  double rs;
  double dx;
  double err_min;
  double a_mix;
  double mu_lo;
  double mu_hi;
  double mu;
  double xmax;
  int nl;
  int nx;
  int nIter;


} input;

void solveSTLS(input in);

void alloc_stls_arrays(input in, double **xx, double **phi,
		       double **AA,  double **GG, double **GG_new,
		       double **SS, double **SSHF);

void free_stls_arrays(double *xx, double *phi,
                      double *AA,  double *GG, double *GG_new,
                      double *SS, double *SSHF);

double compute_mu(input in);

double normalization_condition(double mu, void *pp);

void wave_vector_grid(double *xx, input in);

void compute_phi(double *phi, double *xx, input in);

void compute_phil(double *phil, double *xx, int ll, input in);

double phixl(double yy, void* pp);

double phix0(double yy, void* pp);

int idx2(int xx, int yy, int x_size);

double ssfHF(double yy, void* pp);

void compute_ssfHF(double *SS,  double *xx, input in);

double csch2(double x);

double coth(double x);

void compute_AA(double *AA, double *xx,  input in);

double Axl2(double xx, int ll, input in);

void compute_ssf(double *SS, double *SSHF, double *AA,
                 double *GG, double *phi, double *xx,
                 input in);

void compute_slfc(double *GG, double *SS, 
		  double *xx, input in);

double slfc(double yy, void* pp);

double compute_internal_energy(double *SS, double *xx,
                               input in);

double uex(double yy, void* pp);

void write_text(double *SS, double *GG, 
		double *xx, input in );

void write_bin(double *phi, double *SSHF, double *AA, input in);

void read_text(double *SS, double *GG, 
	       double *xx, input in);

void read_bin(input *in, double **xx, double **phi,
              double **AA, double **GG, double **GG_new,
              double **SS, double **SSHF);

#endif
