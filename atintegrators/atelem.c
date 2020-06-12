/*
 * The file to be included by 'new-style' integrators that support both
 * Matlab and Python.
 */
#ifndef ATELEM_C
#define ATELEM_C

#include "atcommon.h"

/*----------------------------------------------------*/
/*            For the integrator code                 */
/*----------------------------------------------------*/

#if defined(MATLAB_MEX_FILE)

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

static double* atGetDoubleArraySz(const mxArray *ElemData, const char *fieldname, int *msz, int *nsz)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    if (!field) mexErrMsgIdAndTxt("AT:WrongArg", "The required attribute %s is missing.", fieldname);
    *msz = mxGetM(field);  /*Number of rows in the 2-D array*/
    *nsz = mxGetN(field);  /*Number of columns in the 2-D array.*/
    return mxGetDoubles(field);
}

static double* atGetDoubleArray(const mxArray *ElemData, const char *fieldname)
{
    int msz, nsz;
    return atGetDoubleArraySz(ElemData, fieldname, &msz, &nsz);
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

static double* atGetOptionalDoubleArraySz(const mxArray *ElemData, const char *fieldname, int *msz, int *nsz)
{
    double *ptr = NULL;
    mxArray *field=mxGetField(ElemData,0,fieldname);
    if (field) {
        *msz = mxGetM(field);
        *nsz = mxGetN(field);
        ptr = mxGetDoubles(field);
    }
    return ptr;
}

static double* atGetOptionalDoubleArray(const mxArray *ElemData, const char *fieldname)
{
    int msz, nsz;
    return atGetOptionalDoubleArraySz(ElemData, fieldname, &msz, &nsz);
}

#endif /* MATLAB_MEX_FILE */

/*----------------------------------------------------*/
/*            For the Python interface                */
/*----------------------------------------------------*/

#if defined(PYAT)

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

static double *atGetDoubleArraySz(const PyObject *element, char *name, int *msz, int *nsz)
{
    char errmessage[60];
    int ndims;
    npy_intp *dims;
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
    if ((PyArray_FLAGS(array) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
        snprintf(errmessage, 60, "The attribute %s is not Fortran-aligned.", name);
        PyErr_SetString(PyExc_RuntimeError, errmessage);
        return NULL;
    }
    ndims = PyArray_NDIM(array);
    dims = PyArray_SHAPE(array);
    *nsz = (ndims >= 2) ? dims[1] : 0;
    *msz = (ndims >= 1) ? dims[0] : 0;
    return (double *) PyArray_DATA(array);
}

static double *atGetDoubleArray(const PyObject *element, char *name)
{
    int msz, nsz;
    return atGetDoubleArraySz(element, name, &msz, &nsz);
}

static double *atGetOptionalDoubleArraySz(const PyObject *element, char *name, int *msz, int *nsz)
{
    PyObject *obj = PyObject_GetAttrString((PyObject *)element, name);
    if (obj == NULL) {
        PyErr_Clear();
        return NULL;
    }
    return atGetDoubleArraySz(element, name, msz, nsz);
}

static double *atGetOptionalDoubleArray(const PyObject *element, char *name)
{
    int msz, nsz;
    return atGetOptionalDoubleArraySz(element, name, &msz, &nsz);
}

#endif /* defined(PYAT) */

#if defined(PYAT) || defined(MATLAB_MEX_FILE)
#include "attypes.h"

#ifdef __cplusplus
#define C_LINK extern "C"
#else
#define C_LINK
#endif

C_LINK ExportMode struct elem *trackFunction(const atElem *ElemData, struct elem *Elem, double *r_in,
                                      int num_particles, struct parameters *Param);

#endif /* defined(PYAT) || defined(MATLAB_MEX_FILE) */

#endif /*ATELEM_C*/
