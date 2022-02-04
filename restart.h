#ifndef RESTART_H
#define RESTART_H

#include "read_input.h"

// -------------------------------------------------------------------
// FUNCTION USED TO WRITE BINARY FILES FOR GUESS (OR RESTART) STARTING
// FROM TEXT FILES OBTAINED FROM A SIMULATION
// -------------------------------------------------------------------

void create_restart(input in);

// -------------------------------------------------------------------
// FUNCTION USED TO INFER THE FORMAT OF THE TEXT FILES
// -------------------------------------------------------------------

void get_restart_data_format(char * file_name, int *n_lines, int *n_columns);

// -------------------------------------------------------------------
// FUNCTION USED TO SET THE PARAMETERS FOR THE GRID SIZE AND FOR THE
// NUMBER OF MATSUBARA FREQUENCIES
// -------------------------------------------------------------------

void set_nx_nl(int nl1, int nl2, int nc1, int nc2, input *in);

// -------------------------------------------------------------------
// FUNCTION USED TO READ THE TEXT FILES
// -------------------------------------------------------------------

void get_restart_data(char *file_name, int n_lines, int n_columns,
		      double *data, double *xx, input *in);

#endif
