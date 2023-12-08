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

int PyInterface::getOptionalInt(const std::string& name,int defaultValue) {

  try {
    return getInt(name);
  } catch (string&) {
    return defaultValue;
  }

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

AT_FLOAT PyInterface::getOptionalDouble(const std::string& name, AT_FLOAT defaultValue) {

  try {
    return getDouble(name);
  } catch (string&) {
    return defaultValue;
  }

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

  return (AT_FLOAT *) PyArray_DATA(array);

}

AT_FLOAT *PyInterface::getDoubleArray(const std::string& name, std::vector<int64_t> expectedShape) {

  std::vector<int64_t> shape;
  AT_FLOAT *array = getNativeDoubleArray(name,shape);

  if(shape.size()!=expectedShape.size()) {
    Py_DECREF(array);
    throw string(name + " array attribute has shape "+ getShapeStr(shape) +" but " +
           getShapeStr(expectedShape) + " expected");
  }

  bool ok=true;
  size_t d = 0;
  while(ok && d<expectedShape.size()) {
    ok = shape[d] == expectedShape[d];
    d++;
  }
  if( !ok ) {
    Py_DECREF(array);
    throw string(name + " array attribute has shape "+ getShapeStr(shape) +" but " +
           getShapeStr(expectedShape) + " expected");
  }

  return array;

}

AT_FLOAT *PyInterface::getOptionalDoubleArray(const std::string &name, std::vector<int64_t> expectedShape) {

  try {
    return getDoubleArray(name,expectedShape);
  } catch (string&) {
    return nullptr;
  }

}