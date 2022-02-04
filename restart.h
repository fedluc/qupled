#ifndef RESTART_H
#define RESTART_H

#include "read_input.h"

// ---------------------------------------------------------------------
// FUNCTION USED TO CONSTRUCT THE RESTART FILES FROM SIMULATION RESULTS
// ---------------------------------------------------------------------

void create_restart(input in);

void get_restart_data_format(char * file_name, int *n_lines, int *n_columns);

void set_nx_nl(int nl1, int nl2, int nc1, int nc2, input *in);

void get_restart_data(char *file_name, int n_lines, int n_columns,
		      double *data, double *xx, input *in);

#endif
