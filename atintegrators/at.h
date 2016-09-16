#ifndef AT_H
#define AT_H

#ifdef MATLAB_MEX_FILE

#include "mex.h"
#include <matrix.h>

#else

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef NAN
static const double dnan = 0.0 / 0.0;
#define NAN dnan
#endif
#ifndef INFINITY
static const double pinf = 1.0 / 0.0;
#define INFINITY pinf
#endif

#define mxIsFinite isfinite
#define mxIsNaN isnan
#define mxGetNaN() (NAN)
#define mxGetInf() (INFINITY)
#define mxMalloc malloc
#define mxCalloc calloc
#define mxFree free

#if defined __SUNPRO_C
#include <ieeefp.h>
#define isfinite finite
#endif

#ifdef __MWERKS__
#endif

#endif /*MATLAB_MEX_FILE*/

struct parameters
{
  int mode;
  int nturn;
  double RingLength;
  double T0;
};

#endif /*AT_H*/


