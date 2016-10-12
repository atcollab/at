#ifndef ATELEM_C
#define ATELEM_C

#include <stdlib.h>
#include <math.h>
#include "attypes.h"

#if defined(PCWIN) || defined(_WIN32)
#define ExportMode __declspec(dllexport)
#else
#define ExportMode
#endif

/*----------------------------------------------------*/
/*            For the integrator code                 */
/*----------------------------------------------------*/

#if defined(MATLAB_MEX_FILE)

#include <mex.h>
#include <matrix.h>
#define atIsFinite mxIsFinite
#define atIsNaN mxIsNaN
#define atGetNaN mxGetNaN
#define atGetInf mxGetInf
#define atFree mxFree

static void *atMalloc(size_t size)
{
    void *ptr = mxMalloc(size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
}

static void *atCalloc(size_t count, size_t size)
{
    void *ptr = mxCalloc(count, size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
}

#else /* !defined(MATLAB_MEX_FILE) */

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

#define atIsFinite isfinite
#define atIsNaN isnan
#define atGetNaN() (NAN)
#define atGetInf() (INFINITY)
#define atMalloc malloc
#define atCalloc calloc
#define atFree free

#endif /* MATLAB_MEX_FILE */

/*----------------------------------------------------*/
/*            For the Matlab interface                */
/*----------------------------------------------------*/

#if defined(MATLAB_MEX_FILE)

typedef mxArray atElem;
#define check_error()

static long atGetLong(const mxArray *ElemData, const char *fieldname)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    if (!field) mexErrMsgIdAndTxt("AT:WrongArg", "The required attribute %s is missing.", fieldname);
    return (long)mxGetScalar(field);
}

static double atGetDouble(const mxArray *ElemData, const char *fieldname)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    if (!field) mexErrMsgIdAndTxt("AT:WrongArg", "The required attribute %s is missing.", fieldname);
    return mxGetScalar(field);
}

static double* atGetDoubleArray(const mxArray *ElemData, const char *fieldname)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    if (!field) mexErrMsgIdAndTxt("AT:WrongArg", "The required attribute %s is missing.", fieldname);
    return mxGetPr(field);
}

static long atGetOptionalLong(const mxArray *ElemData, const char *fieldname, long default_value)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    return (field) ? (long)mxGetScalar(field) : default_value;
}

static double atGetOptionalDouble(const mxArray *ElemData, const char *fieldname, double default_value)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    return (field) ? mxGetScalar(field) : default_value;
}

static double* atGetOptionalDoubleArray(const mxArray *ElemData, const char *fieldname)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    return (field) ? mxGetPr(field) : NULL;
}

#endif /* MATLAB_MEX_FILE */

/*----------------------------------------------------*/
/*            For the Python interface                */
/*----------------------------------------------------*/

#if defined(PYAT)

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#if PY_MAJOR_VERSION >= 3
#define NUMPY_IMPORT_ARRAY_RETVAL NULL
#define NUMPY_IMPORT_ARRAY_TYPE void *
#else
#define NUMPY_IMPORT_ARRAY_RETVAL
#define NUMPY_IMPORT_ARRAY_TYPE void
#define PyLong_AsLong PyInt_AsLong
#endif

typedef PyObject atElem;
#define check_error() if (PyErr_Occurred()) return NULL

static int array_imported = 0;

static NUMPY_IMPORT_ARRAY_TYPE init_numpy(void)
{
    import_array();
    return NUMPY_IMPORT_ARRAY_RETVAL;
}

static long atGetLong(const PyObject *element, const char *name)
{
    const PyObject *attr = PyObject_GetAttrString((PyObject *)element, name);
    if (!attr) return 0L;
    return PyLong_AsLong((PyObject *)attr);
}

static double atGetDouble(const PyObject *element, const char *name)
{
    const PyObject *attr = PyObject_GetAttrString((PyObject *)element, name);
    if (!attr) return 0.0;
    return PyFloat_AsDouble((PyObject *)attr);
}

static long atGetOptionalLong(const PyObject *element, const char *name, long default_value)
{
    long l = atGetLong(element, name);
    if (PyErr_Occurred()) {
        PyErr_Clear();
        l = default_value;
    }
    return l;
}

static double atGetOptionalDouble(const PyObject *element, const char *name, double default_value)
{
    double d = atGetDouble(element, name);
    if (PyErr_Occurred()) {
        PyErr_Clear();
        d = default_value;
    }
    return d;
}

static double *atGetDoubleArray(const PyObject *element, char *name)
{
    char errmessage[60];
    PyArrayObject *array;
    if (!array_imported) {
        init_numpy();
        array_imported = 1;
    }
    array = (PyArrayObject *) PyObject_GetAttrString((PyObject *)element, name);
    if (array == NULL) {
        return NULL;
    }
    if (!PyArray_Check(array)) {
        snprintf(errmessage, 60, "The attribute %s is not an array.", name);
        PyErr_SetString(PyExc_RuntimeError, errmessage);
        return NULL;
    }
    if (PyArray_TYPE(array) != NPY_DOUBLE) {
        snprintf(errmessage, 60, "The attribute %s is not a double array.", name);
        PyErr_SetString(PyExc_RuntimeError, errmessage);
        return NULL;
    }
    if ((PyArray_FLAGS(array) & NPY_ARRAY_CARRAY_RO) != NPY_ARRAY_CARRAY_RO) {
        snprintf(errmessage, 60, "The attribute %s is not C-aligned.", name);
        PyErr_SetString(PyExc_RuntimeError, errmessage);
        return NULL;
    }
    return PyArray_DATA(array);
}

static double *atGetOptionalDoubleArray(const PyObject *element, char *name)
{
    PyObject *obj = PyObject_GetAttrString((PyObject *)element, name);
    if (obj == NULL) {
        PyErr_Clear();
        return NULL;
    }
    return atGetDoubleArray(element, name);
}

#endif /* defined(PYAT) */

ExportMode struct elem *trackFunction(const atElem *ElemData, struct elem *Elem, double *r_in,
                                      int num_particles, struct parameters *Param);

#endif /*ATELEM_C*/
