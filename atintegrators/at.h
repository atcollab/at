#ifndef AT_H
#define AT_H

#include <stdlib.h>
#include <math.h>
#include "attypes.h"

#ifdef MATLAB_MEX_FILE

#include "mex.h"
#include <matrix.h>

#else

#if defined(_WIN32)

#include <Windows.h>
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#define isfinite(x) _finite(x)
/* See https://blogs.msdn.microsoft.com/oldnewthing/20100305-00/?p=14713 */
DECLSPEC_SELECTANY extern const float FLOAT_NaN = ((float)((1e308 * 10)*0.));
#define NAN FLOAT_NaN
DECLSPEC_SELECTANY extern const float FLOAT_POSITIVE_INFINITY = ((float)(1e308 * 10));
#define INFINITY FLOAT_POSITIVE_INFINITY
typedef int bool;
#define false 0
#define true 1

#else /* !defined(_WIN32) */

#if defined __SUNPRO_C
#include <ieeefp.h>
#define isfinite finite
#endif
#ifndef NAN
static const double dnan = 0.0 / 0.0;
#define NAN dnan
#endif
#ifndef INFINITY
static const double pinf = 1.0 / 0.0;
#define INFINITY pinf
#endif
#include <stdbool.h>

#endif /* defined(_WIN32) */

#define mxIsFinite isfinite
#define mxIsNaN isnan
#define mxGetNaN() (NAN)
#define mxGetInf() (INFINITY)
#define mxMalloc malloc
#define mxCalloc calloc
#define mxFree free

#endif /*MATLAB_MEX_FILE*/

#endif /*AT_H*/


