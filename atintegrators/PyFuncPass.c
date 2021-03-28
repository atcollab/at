#define PY_SSIZE_T_CLEAN
#include <Python.h> 
#include "atelem.c"
#include "numpy/arrayobject.h"

struct elem
{
    double Length;
    PyObject *function;
    PyObject *kwargs;
};

static PyObject *GetFunction(const atElem *ElemData)
{
  PyObject *eName = PyObject_GetAttrString((PyObject *)ElemData, "FamName");
  PyObject *pName = PyObject_GetAttrString((PyObject *)ElemData, "pyModule");
  PyObject *fName = PyObject_GetAttrString((PyObject *)ElemData, "pyFun");
  const char *s1 = PyUnicode_AsUTF8(eName);
  const char *s2 = PyUnicode_AsUTF8(pName);
  const char *s3 = PyUnicode_AsUTF8(fName);
  PyObject *pModule = PyImport_Import(pName);
  if (!pModule){
      printf("PyFuncPass: could not import pyModule %s in %s \n", s2, s1);
    }
  PyObject *function = PyObject_GetAttrString(pModule, PyUnicode_AsUTF8(fName));
  if (!function){
      printf("PyFuncPass: could not import pyFun %s in %s \n", s3, s1);
    }
  if (!PyCallable_Check(function)){
      printf("PyFuncPass: pyFun %s in %s not callable \n", s3, s1);
    }
  return function;
}


static PyObject *Buildkwargs(const atElem *ElemData)
{
  PyObject *kwargs = PyDict_New();
  PyDict_SetItemString(kwargs,(char *)"elem",(PyObject *)ElemData);
  return kwargs;  
}


static PyObject *Buildargs(double *r_in, int num_particles)
{
  npy_intp outdims[1];
  outdims[0] = 6*num_particles;
  PyObject *rin = PyArray_SimpleNewFromData(1, outdims, NPY_DOUBLE, r_in);
  if (!rin){
      printf("PyFuncPass: could not generate pyArray rin");
    }
  return PyTuple_Pack(1,rin);
}


void PyFuncPass(double *r_in, PyObject *function, PyObject *kwargs, int num_particles)
{   
  PyObject *args = Buildargs(r_in, num_particles);
  PyObject_Call(function, args, kwargs);
}


#if defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        _import_array();
        double Length;
        Length=atGetDouble(ElemData,"Length"); check_error();  
        PyObject *function = GetFunction(ElemData); check_error();
        PyObject *kwargs = Buildkwargs(ElemData); check_error();      
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->function=function;
        Elem->kwargs=kwargs;
    } 
    PyFuncPass(r_in, Elem->function, Elem->kwargs, num_particles);
    return Elem;
}

MODULE_DEF(PyFuncPass)        /* Dummy module initialisation */

#endif
