#include "PyInterface.h"
#include <string>

using namespace std;

void PyInterface::setObject(PyObject *obj) {
  self = obj;
}

std::string PyInterface::getShapeStr(size_t nDim, int64_t *sizes) {
  string ret = "(";
  for(size_t i=nDim;i<nDim;i++) {
    ret.append(to_string(sizes[i]));
    if(i<nDim-1) ret.append(",");
  }
  ret.append(")");
  return ret;
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
    throw(name + " attribute not found");
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

AT_FLOAT *PyInterface::getDoubleArray(const std::string& name, std::vector<int64_t> expectedShape) {

  PyArrayObject *array = (PyArrayObject *) PyObject_GetAttrString(self, name.c_str());
  if (array == nullptr)
    throw(name + " array attribute not found");

  if ((PyArray_FLAGS(array) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
    Py_DECREF(array);
    throw(name + " array attribute is not Fortran-aligned");
  }

  if (PyArray_TYPE(array) != NPY_DOUBLE) {
    Py_DECREF(array);
    throw(name + " array attribute is not a double array");
  }

  size_t nDim = PyArray_NDIM(array);
  int64_t *dims = PyArray_SHAPE(array);
  if(nDim!=expectedShape.size()) {
    Py_DECREF(array);
    throw (name + " array attribute has shape "+ getShapeStr(nDim,dims) +" but " +
           getShapeStr(expectedShape.size(),expectedShape.data()) + " expected");
  }

  bool ok=true;
  size_t d = 0;
  while(ok && d<expectedShape.size()) {
    ok = dims[d] == expectedShape[d];
    d++;
  }
  if( !ok ) {
    Py_DECREF(array);
    throw (name + " array attribute has shape "+ getShapeStr(nDim,dims) +" but " +
           getShapeStr(expectedShape.size(),expectedShape.data()) + " expected");
  }

  return (double *) PyArray_DATA(array);

}

AT_FLOAT *PyInterface::getOptionalDoubleArray(const std::string &name, std::vector<int64_t> expectedShape) {

  try {
    return getDoubleArray(name,expectedShape);
  } catch (string&) {
    return nullptr;
  }

}