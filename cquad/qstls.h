#ifndef QSTLS_H
#define QSTLS_H

#include "solvers.h"

// -------------------------------------------------------------------
// FUNCTION USED TO COMPUTE THE ...
// -------------------------------------------------------------------

void compute_psi_xlw(double *psi_xlw, double *xx, input in);

double psi_x0w_t(double uu, void* pp);
double psi_x0w_q(double uu, void* pp);
double psi_xlw_t(double uu, void* pp);
double psi_xlw_q(double uu, void* pp);


void compute_ssf_qstls(double *SS, double *SSHF, double *psi,
                       double *phi, double *xx, input in);

void compute_psi(double *psi, double *psi_xlw, double *SS,
		 double *xx, input in);

double psiw(double ww, void* pp);

int idx3(int xx, int yy, int zz,
         int x_size, int y_size);

#endif
