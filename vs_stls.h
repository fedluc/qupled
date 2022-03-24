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
  double *el[VST_NUMEL];
  
} vs_thermo;

#endif

