#ifndef ATCOMMON_H
#define ATCOMMON_H

/* All builds */
#include <stdlib.h>
#include <math.h>
#include "attypes.h"

/* All Windows builds */
#if defined(PCWIN) || defined(_WIN32)
#define ExportMode __declspec(dllexport)
#else
#define ExportMode
#endif


#ifdef MATLAB_MEX_FILE
/* Matlab only */
#include <mex.h>
#include <matrix.h>

#else

#if defined(_WIN32) && (_MSC_VER < 1800)
/* Python Windows builds */
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
#include <stdbool.h>
#if defined __SUNPRO_C
/* Python Sun builds? */
#include <ieeefp.h>
#define isfinite finite
#endif
/* All other Python builds */
#ifndef NAN
static const double dnan = 0.0 / 0.0;
#define NAN dnan
#endif
#ifndef INFINITY
static const double pinf = 1.0 / 0.0;
#define INFINITY pinf
#endif

#endif /* defined(_WIN32) */

#endif /*MATLAB_MEX_FILE*/

#endif /*ATCOMMON_H*/
