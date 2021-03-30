#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>
#include "atelem.c"

#define NUMPY_IMPORT_ARRAY_RETVAL NULL
#define NUMPY_IMPORT_ARRAY_TYPE void *

struct elem
{
    double Length;
    PyObject *pyfunction;
    PyObject *kwargs;
};

static PyObject *GetFunction(const atElem *ElemData)
{
  PyObject *eName;
  PyObject *pName;
  PyObject *fName;
  PyObject *pModule;
  PyObject *pyfunction;
  const char *s1;
  const char *s2;
  const char *s3;

  eName = PyObject_GetAttrString((PyObject *)ElemData, "FamName");
  pName = PyObject_GetAttrString((PyObject *)ElemData, "pyModule");
  fName = PyObject_GetAttrString((PyObject *)ElemData, "pyFun");
  s1 = PyUnicode_AsUTF8(eName);
  s2 = PyUnicode_AsUTF8(pName);
  s3 = PyUnicode_AsUTF8(fName);
  pModule = PyImport_Import(pName);
  if (!pModule){
      printf("PyFuncPass: could not import pyModule %s in %s \n", s2, s1);
    }
  pyfunction = PyObject_GetAttrString(pModule, PyUnicode_AsUTF8(fName));
  if (!pyfunction){
      printf("PyFuncPass: could not import pyFun %s in %s \n", s3, s1);
    }
  if (!PyCallable_Check(pyfunction)){
      printf("PyFuncPass: pyFun %s in %s not callable \n", s3, s1);
    }
  Py_DECREF(eName);
  Py_DECREF(pName);
  Py_DECREF(fName);
  Py_DECREF(pModule);
  return pyfunction;
}


static PyObject *Buildkwargs(const atElem *ElemData)
{
  PyObject *kwargs;
  kwargs = PyDict_New();
  PyDict_SetItemString(kwargs,(char *)"elem",(PyObject *)ElemData);
  return kwargs;  
}


static PyObject *Buildargs(double *r_in, int num_particles)
{
  npy_intp outdims[1];
  outdims[0] = 6*num_particles;
  PyObject *rin;
  rin = PyArray_SimpleNewFromData(1, outdims, NPY_DOUBLE, r_in);
  if (!rin){
      printf("PyFuncPass: could not generate pyArray rin");
    }
  return PyTuple_Pack(1,rin);
}


void PyFuncPass(double *r_in, PyObject *function, PyObject *kwargs, int num_particles)
{   
  PyObject *args;
  args = Buildargs(r_in, num_particles);
  PyObject_Call(function, args, kwargs);
  Py_DECREF(args);
}


#if defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        _import_array();
        double Length;
        PyObject *pyfunction;
        PyObject *kwargs;
        Length=atGetDouble(ElemData,"Length"); check_error();  
        pyfunction = GetFunction(ElemData);
        kwargs = Buildkwargs(ElemData);      
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->pyfunction=pyfunction;
        Elem->kwargs=kwargs;
    } 
    PyFuncPass(r_in, Elem->pyfunction, Elem->kwargs, num_particles);
    return Elem;
}

MODULE_DEF(PyFuncPass)        /* Dummy module initialisation */

#endif
