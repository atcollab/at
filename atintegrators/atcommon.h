#ifndef ATCOMMON_H
#define ATCOMMON_H

#ifdef PYAT
/* Python.h must be included first. */
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#define STR(name) XSTR(name)
#define XSTR(name) #name
/* Handle differences between Python 2 and Python 3 in declaring extension
   modules and the numpy API. */
#if PY_MAJOR_VERSION >= 3
#define MOD_INIT_STR(name) PyInit_##name(void)
#define MOD_ERROR_VAL NULL
#define MOD_SUCCESS_VAL(val) val
#define NUMPY_IMPORT_ARRAY_RETVAL NULL
#define NUMPY_IMPORT_ARRAY_TYPE void *
#else
#define MOD_INIT_STR(name) init##name(void)
#define MOD_ERROR_VAL
#define MOD_SUCCESS_VAL(val)
#define NUMPY_IMPORT_ARRAY_RETVAL
#define NUMPY_IMPORT_ARRAY_TYPE void
#define PyLong_AsLong PyInt_AsLong
#endif /*PY_MAJOR_VERSION*/
#define MOD_INIT(name) MOD_INIT_STR(name)


#if defined(_WIN32)    /* Create a dummy module initialisation function for Windows */
#define MODULE_DEF(name) PyMODINIT_FUNC MOD_INIT(name) {return MOD_ERROR_VAL;}
#endif /* _WIN32 */
#endif /* PYAT */

#ifndef MODULE_DEF
#define MODULE_DEF(name)
#endif

/* All builds */
#include <stdlib.h>
#include <math.h>

/* All Windows builds */
#if defined(PCWIN) || defined(_WIN32)
#define ExportMode __declspec(dllexport)
#else
#define ExportMode
#endif


#ifdef MATLAB_MEX_FILE
/* Matlab only */
#include <mex.h>
#ifndef OCTAVE
#include <matrix.h>
#endif

/* Get ready for R2018a C matrix API */
#ifndef mxGetDoubles
#define mxGetDoubles mxGetPr
typedef double mxDouble;
#endif

#else

#if defined(_WIN32) && (_MSC_VER < 1800)
/* Python Windows builds */
#include <Windows.h>
#include <float.h>
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#define isfinite(x) _finite(x)
/* See https://blogs.msdn.microsoft.com/oldnewthing/20100305-00/?p=14713 */
DECLSPEC_SELECTANY extern const float FLOAT_NaN = ((float)((1e308 * 10)*0.));
#define NAN FLOAT_NaN
DECLSPEC_SELECTANY extern const float FLOAT_POSITIVE_INFINITY = ((float)(1e308 * 10));
#define INFINITY FLOAT_POSITIVE_INFINITY
#ifndef __cplusplus
typedef int bool;
#endif
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
