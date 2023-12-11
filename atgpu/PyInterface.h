#ifndef AT_GPU_PYINTERFACE_H
#define AT_GPU_PYINTERFACE_H
#include "AbstractInterface.h"
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

class PyInterface: public AbstractInterface {

public:

  std::string getString(const std::string& name);
  int getInt(const std::string& name);
  AT_FLOAT getDouble(const std::string& name);
  AT_FLOAT *getDoubleArray(const std::string& name,std::vector<int64_t> expectedShape);
  AT_FLOAT *getNativeDoubleArray(const std::string& name,std::vector<int64_t>& shape);

  void setObject(PyObject *obj);

private:

  PyObject *self;

};

#endif //AT_GPU_PYINTERFACE_H
