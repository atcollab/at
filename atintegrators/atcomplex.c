#ifndef ATCOMPLEX_C
#define ATCOMPLEX_C

#ifndef __cplusplus
/* Complex attributes for C mex-files */
#include <complex.h>

#if defined(MATLAB_MEX_FILE)
/*----------------------------------------------------*/
/*            For the Matlab interface                */
/*----------------------------------------------------*/

static double _Complex atGetComplex(const mxArray *ElemData, const char *fieldname)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    if (!field) mexErrMsgIdAndTxt("AT:WrongArg", "The required attribute %s is missing.", fieldname);
    return *(double _Complex *) mxGetComplexDoubles(field);
}

static double _Complex* atGetComplexArraySz(const mxArray *ElemData, const char *fieldname, int *msz, int *nsz)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    if (!field) mexErrMsgIdAndTxt("AT:WrongArg", "The required attribute %s is missing.", fieldname);
    *msz = mxGetM(field);  /*Number of rows in the 2-D array*/
    *nsz = mxGetN(field);  /*Number of columns in the 2-D array.*/
    return (double _Complex *)mxGetComplexDoubles(field);
}

static double _Complex* atGetComplexArray(const mxArray *ElemData, const char *fieldname)
{
    int msz, nsz;
    return atGetComplexArraySz(ElemData, fieldname, &msz, &nsz);
}

static double _Complex atGetOptionalComplex(const mxArray *ElemData, const char *fieldname, double _Complex default_value)
{
    mxArray *field=mxGetField(ElemData,0,fieldname);
    return (field) ? *(double _Complex *) mxGetComplexDoubles(field) : default_value;
}

static double _Complex* atGetOptionalComplexArraySz(const mxArray *ElemData, const char *fieldname, int *msz, int *nsz)
{
    mxComplexDouble *ptr = NULL;
    mxArray *field=mxGetField(ElemData,0,fieldname);
    if (field) {
        *msz = mxGetM(field);
        *nsz = mxGetN(field);
        ptr = mxGetComplexDoubles(field);
    }
    return (double _Complex *) ptr;
}

static double _Complex* atGetOptionalComplexArray(const mxArray *ElemData, const char *fieldname)
{
    int msz, nsz;
    return atGetOptionalComplexArraySz(ElemData, fieldname, &msz, &nsz);
}

#endif /* MATLAB_MEX_FILE */

#if defined(PYAT)
/*----------------------------------------------------*/
/*            For the Python interface                */
/*----------------------------------------------------*/

static double _Complex *atGetComplexArraySz(const PyObject *element, char *name, int *msz, int *nsz)
{
    return (double _Complex *) atGetAnyArraySz(element, name, NPY_CDOUBLE, "complex", msz, nsz);
}

static double _Complex *atGetComplexArray(const PyObject *element, char *name)
{
    int msz, nsz;
    return atGetComplexArraySz(element, name, &msz, &nsz);
}

static double _Complex *atGetOptionalComplexArraySz(const PyObject *element, char *name, int *msz, int *nsz)
{
    PyObject *obj = PyObject_GetAttrString((PyObject *)element, name);
    if (obj == NULL) {
        PyErr_Clear();
        return NULL;
    }
    return atGetComplexArraySz(element, name, msz, nsz);
}

static double _Complex *atGetOptionalComplexArray(const PyObject *element, char *name)
{
    int msz, nsz;
    return atGetOptionalComplexArraySz(element, name, &msz, &nsz);
}

#endif /* defined(PYAT) */
#endif /*__cplusplus*/

#endif /*ATCOMPLEX_C*/
