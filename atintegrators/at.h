#ifndef AT_H
#define AT_H

#ifdef MATLAB_MEX_FILE

#include <matrix.h>

#else

#include <math.h>
#include <stdlib.h>

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
#define mxFree free

#if defined __SUNPRO_C
#include <ieeefp.h>
#define isfinite finite
#endif

#ifdef __MWERKS__
#endif

#endif /*MATLAB_MEX_FILE*/

#endif /*AT_H*/
