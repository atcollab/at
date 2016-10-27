/*
 * The header file used when defining an 'old-style' Matlab integrator.
 * 'New-style' integrators should support both Matlab and Python and
 * should include "atelem.c".
 */
#ifndef AT_H
#define AT_H

#include "atcommon.h"

#ifndef MATLAB_MEX_FILE
#define mxIsFinite isfinite
#define mxIsNaN isnan
#define mxGetNaN() (NAN)
#define mxGetInf() (INFINITY)
#define mxMalloc malloc
#define mxCalloc calloc
#define mxFree free
#endif

#endif /*AT_H*/
