/*
 * This file provides function returning configuration variables
 * of the C compiler
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject *isopenmp(PyObject *self)
{
#ifdef _OPENMP
    Py_RETURN_TRUE;
#else
    Py_RETURN_FALSE;
#endif /*_OPENMP)*/
}


static PyObject *ismpi(PyObject *self)
{
#ifdef MPI
    Py_RETURN_TRUE;
#else
    Py_RETURN_FALSE;
#endif /*MPI)*/
}

static PyObject *iscuda(PyObject *self)
{
#ifdef CUDA
  Py_RETURN_TRUE;
#else
  Py_RETURN_FALSE;
#endif /*CUDA)*/
}

static PyObject *isopencl(PyObject *self)
{
#ifdef OPENCL
  Py_RETURN_TRUE;
#else
  Py_RETURN_FALSE;
#endif /*OPENCL)*/
}

static PyMethodDef methods[] = {
    {"isopenmp",  (PyCFunction)isopenmp, METH_NOARGS,
    PyDoc_STR("isopenmp()\n\n"
              "Return whether OpenMP is active.\n"
             )},
    {"ismpi",  (PyCFunction)ismpi, METH_NOARGS,
    PyDoc_STR("ismpi()\n\n"
              "Return whether MPI is active.\n"
             )},
    {"iscuda",  (PyCFunction)iscuda, METH_NOARGS,
            PyDoc_STR("iscuda()\n\n"
                      "Return whether CUDA is active.\n"
            )},
    {"isopencl",  (PyCFunction)isopencl, METH_NOARGS,
            PyDoc_STR("isopencl()\n\n"
                      "Return whether OpenCL is active.\n"
            )},
   {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC PyInit_cconfig(void)
{
    static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "cconfig",    /* m_name */
    PyDoc_STR("C configuration"),      /* m_doc */
    -1,           /* m_size */
    methods,    /* m_methods */
    NULL,         /* m_slots */
    NULL,         /* m_traverse */
    NULL,         /* m_clear */
    NULL,         /* m_free */
    };
    return PyModule_Create(&moduledef);
}
