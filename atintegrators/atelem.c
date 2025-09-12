/*
 * The file to be included by 'new-style' integrators that support both
 * Matlab and Python.
 */
#ifndef ATELEM_C
#define ATELEM_C

#include "atcommon.h"
#include "attypes.h"
#include "atconstants.h"

/*----------------------------------------------------*/
/*            For the integrator code                 */
/*----------------------------------------------------*/

#if defined(MATLAB_MEX_FILE)

#include <string.h>
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

#define SQR(X) ((X)*(X))

/* AT coordinates */

#ifndef __cplusplus
#define x_ 0
#define px_ 1
#define y_ 2
#define py_ 3
#define delta_ 4
#define ct_ 5
#endif /*__cplusplus */

/*----------------------------------------------------*/
/*            For the Matlab interface                */
/*----------------------------------------------------*/

#if defined(MATLAB_MEX_FILE)

typedef mxArray atElem;
#define check_error()
#define atError(...) mexErrMsgIdAndTxt("AT:PassError", __VA_ARGS__)
#define atWarning(...) mexWarnMsgIdAndTxt("AT:PassWarning", __VA_ARGS__)
#define atPrintf(...) mexPrintf(__VA_ARGS__)
#include "ringproperties.c"

double atEnergy(double ringenergy, double elemenergy)
{
    if (ringenergy!=0.0)
        return ringenergy;
    else
        if (elemenergy!=0.0)
            return elemenergy;
        else {
            atError("Energy not defined.");
            return 0.0;   /* Never reached but makes the compiler happy */
        }
}

double atGamma(double ringenergy, double elemenergy, double rest_energy)
{
    double energy = atEnergy(ringenergy, elemenergy);
    if (rest_energy == 0.0)
        return 1.0E-9 * energy / __E0;
    else
        return energy / rest_energy;
}

static mxArray *get_field(const mxArray *pm, const char *fieldname)
{
   mxArray *field;
   if (fieldname[0] == '_') {   /* replace leading '-' by trailing '_' */
      size_t n=strlen(fieldname);
      char *buffer=strcpy((char *)malloc(n+1),fieldname+1);
      buffer[n-1]='_';
      buffer[n]='\0';
      field = mxGetField(pm,0,buffer);
      free(buffer);
   }
   else {
      field = mxGetField(pm,0,fieldname);
   }
   return field;
}

static long atGetLong(const mxArray *ElemData, const char *fieldname)
{
    mxArray *field=get_field(ElemData,fieldname);
    if (!field) mexErrMsgIdAndTxt("AT:WrongArg", "The required attribute %s is missing.", fieldname);
    return (long)mxGetScalar(field);
}

static double atGetDouble(const mxArray *ElemData, const char *fieldname)
{
    mxArray *field=get_field(ElemData,fieldname);
    if (!field) mexErrMsgIdAndTxt("AT:WrongArg", "The required attribute %s is missing.", fieldname);
    return mxGetScalar(field);
}

static double* atGetDoubleArraySz(const mxArray *ElemData, const char *fieldname, int *msz, int *nsz)
{
     mxArray *field=get_field(ElemData,fieldname);
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
    mxArray *field=get_field(ElemData,fieldname);
    return (field) ? (long)mxGetScalar(field) : default_value;
}

static double atGetOptionalDouble(const mxArray *ElemData, const char *fieldname, double default_value)
{
    mxArray *field=get_field(ElemData,fieldname);
    return (field) ? mxGetScalar(field) : default_value;
}

void atCheckArrayDims(const mxArray *ElemData, char *fieldname, int ndim, int *dims)
{
    const mwSize *dptr, *dlim;
    int i;
    mwSize nd;

    mxArray *field=get_field(ElemData,fieldname);
    if (!field)
        mexErrMsgIdAndTxt("AT:WrongArg", "The required attribute %s is missing.", fieldname);
    nd=mxGetNumberOfDimensions(field);
    if (nd != ndim)
        mexErrMsgIdAndTxt("AT:WrongArg", "%s should have %d dimensions instead of %d.", fieldname, ndim, nd);
    dptr = mxGetDimensions(field);
    for (i=0; i < ndim; i++){
        if (dptr[i] != dims[i]){
            mexErrMsgIdAndTxt("AT:WrongArg", "%s dimension %d has size %d instead of %d", fieldname, i, dptr[i], dims[i]);
        }
    }
}

static double* atGetOptionalDoubleArraySz(const mxArray *ElemData, const char *fieldname, int *msz, int *nsz)
{
    double *ptr = NULL;
    mxArray *field=get_field(ElemData,fieldname);
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
#define atError(...) return (struct elem *) PyErr_Format(PyExc_ValueError, __VA_ARGS__)
#define atWarning(...) if (PyErr_WarnFormat(PyExc_RuntimeWarning, 0, __VA_ARGS__) != 0) return NULL
#define atPrintf(...) PySys_WriteStdout(__VA_ARGS__)
#define atEnergy(ringenergy,elemenergy) (ringenergy)
#define atGamma(ringenergy,elemenergy,rest_energy) ((rest_energy) == 0.0 ? 1.0E-9*(ringenergy)/__E0 : (ringenergy)/(rest_energy))

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
    long l = PyLong_AsLong((PyObject *)attr);
    Py_DECREF(attr);
    return l;
}

static double atGetDouble(const PyObject *element, const char *name)
{
    const PyObject *attr = PyObject_GetAttrString((PyObject *)element, name);
    if (!attr) return 0.0;
    double d = PyFloat_AsDouble((PyObject *)attr);
    Py_DECREF(attr);
    return d;
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

void atCheckArrayDims(const PyObject *element, char *name, int ndim, int *dims)
{
    char errmessage[60];
    PyArrayObject *array = (PyArrayObject *) PyObject_GetAttrString((PyObject *)element, name);
    if (!array_imported) {
        init_numpy();
        array_imported = 1;
    }
    Py_DECREF(array);
    if (array == NULL) {
        snprintf(errmessage, 60, "The required attribute %s is missing.", name);
        PyErr_SetString(PyExc_RuntimeError, errmessage);
    }
    int ndima, i;
    npy_intp *dimsa;
    ndima = PyArray_NDIM(array);
    dimsa = PyArray_SHAPE(array);

    if (ndima != ndim) {
        snprintf(errmessage, 60, "The attribute %s should have %d dimensions instead of %d.", name, ndim, ndima);
        PyErr_SetString(PyExc_RuntimeError, errmessage);
    }
    for(i=0; i<ndim;i++){
        if (dims[i] != dimsa[i]) {
            snprintf(errmessage, 60, "The attribute %s dimension %d has the wrong size", name, i);
            PyErr_SetString(PyExc_RuntimeError, errmessage);
        }
    }
}

static double *atGetArrayData(PyArrayObject *array, char *name, int atype, int *msz, int *nsz)
{
    char errmessage[60];
    int ndims;
    npy_intp *dims;
    if (!array_imported) {
        init_numpy();
        array_imported = 1;
    }
    Py_DECREF(array);
    if (!PyArray_Check(array)) {
        snprintf(errmessage, 60, "The attribute %s is not an array.", name);
        PyErr_SetString(PyExc_RuntimeError, errmessage);
        return NULL;
    }
    if (PyArray_TYPE(array) != atype) {
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
    *nsz = (ndims >= 2) ? (int)dims[1] : 0;
    *msz = (ndims >= 1) ? (int)dims[0] : 0;
    return (double *) PyArray_DATA(array);
}

static double *atGetDoubleArraySz(const PyObject *element, char *name, int *msz, int *nsz)
{
    PyArrayObject *array = (PyArrayObject *) PyObject_GetAttrString((PyObject *)element, name);
    if (array == NULL) {
        *msz=0;
        *nsz=0;
        return NULL;
    }
    return (double *) atGetArrayData(array, name, NPY_DOUBLE, msz, nsz);
}

static double *atGetDoubleArray(const PyObject *element, char *name)
{
    int msz, nsz;
    return atGetDoubleArraySz(element, name, &msz, &nsz);
}

static double *atGetOptionalDoubleArraySz(const PyObject *element, char *name, int *msz, int *nsz)
{
    PyArrayObject *array = (PyArrayObject *) PyObject_GetAttrString((PyObject *)element, name);
    if (array == NULL) {
        PyErr_Clear();
        *msz=0;
        *nsz=0;
        return NULL;
    }
    return (double *) atGetArrayData(array, name, NPY_DOUBLE, msz, nsz);
}

static double *atGetOptionalDoubleArray(const PyObject *element, char *name)
{
    int msz, nsz;
    return atGetOptionalDoubleArraySz(element, name, &msz, &nsz);
}

#endif /* defined(PYAT) */

#if defined(PYAT) || defined(MATLAB_MEX_FILE)

#ifdef __cplusplus
#define C_LINK extern "C"
#else
#define C_LINK
#endif

C_LINK ExportMode struct elem *trackFunction(const atElem *ElemData, struct elem *Elem, double *r_in,
                                      int num_particles, struct parameters *Param);

#endif /* defined(PYAT) || defined(MATLAB_MEX_FILE) */

#endif /*ATELEM_C*/
