#ifndef AT_GPU_PYINTERFACE_H
#define AT_GPU_PYINTERFACE_H
#include "AbstractInterface.h"
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

class PyInterface: public AbstractInterface {

public:

  std::string getString(const std::string& name);
  AT_FLOAT getDouble(const std::string& name);
  AT_FLOAT getOptionalDouble(const std::string& name,AT_FLOAT defaultValue);
  AT_FLOAT *getDoubleArray(const std::string& name,std::vector<int64_t> expectedShape);
  AT_FLOAT *getOptionalDoubleArray(const std::string& name,std::vector<int64_t> expectedShape);

  void setObject(PyObject *obj);

private:

  PyObject *self;
  std::string getShapeStr(size_t nDim,int64_t *sizes);

};

#endif //AT_GPU_PYINTERFACE_H
