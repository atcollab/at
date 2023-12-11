#include "PyInterface.h"
#include <string>

using namespace std;

void PyInterface::setObject(PyObject *obj) {
  self = obj;
}

int PyInterface::getInt(const std::string& name) {

  PyObject *attr = PyObject_GetAttrString(self, name.c_str());
  if (!attr)
    throw(name + " attribute not found");
  Py_DECREF(attr);
  return (int)PyLong_AsLong(attr);

}

std::string PyInterface::getString(const std::string& name) {

  PyObject *attr = PyObject_GetAttrString(self, name.c_str());
  if (!attr)
    throw(name + " attribute not found");
  Py_DECREF(attr);
  return PyUnicode_AsUTF8(attr);

}


AT_FLOAT PyInterface::getDouble(const std::string& name) {

  PyObject *attr = PyObject_GetAttrString(self, name.c_str());
  if (!attr)
    throw string(name + " attribute not found");
  Py_DECREF(attr);
  return PyFloat_AsDouble(attr);

}


AT_FLOAT *PyInterface::getNativeDoubleArray(const std::string& name,std::vector<int64_t>& shape) {

  PyArrayObject *array = (PyArrayObject *) PyObject_GetAttrString(self, name.c_str());
  if (array == nullptr)
    throw string(name + " array attribute not found");

  if ((PyArray_FLAGS(array) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
    Py_DECREF(array);
    throw string(name + " array attribute is not Fortran-aligned");
  }

  if (PyArray_TYPE(array) != NPY_DOUBLE) {
    Py_DECREF(array);
    throw string(name + " array attribute is not a double array");
  }

  size_t nDim = PyArray_NDIM(array);
  int64_t *dims = PyArray_SHAPE(array);

  shape.resize(nDim);
  for(int i=0;i<nDim;i++)
    shape[i] = dims[i];

  AT_FLOAT *ptr = (AT_FLOAT *)PyArray_DATA(array);
  Py_DECREF(array);
  return ptr;

}


