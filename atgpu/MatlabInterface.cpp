#include "MatlabInterface.h"
#include "AbstractGPU.h"

using namespace std;

void MatlabInterface::setObject(mxArray *obj) {
  elem = obj;
}

mxArray *MatlabInterface::getField(const mxArray *pm, const std::string& name) {

  mxArray *field;
  if (name[0] == '_') {
    // replace leading '_' by trailing '_'
    string newName = name.substr(1) + '_';
    field = mxGetField(pm,0,newName.c_str());
  } else {
    field = mxGetField(pm,0,name.c_str());
  }
  return field;

}

int MatlabInterface::getInt(const std::string& name) {

  mxArray *field=getField(elem,name);
  if (!field)
    throw string("The required attribute " + name + " is missing.");
  return (int)mxGetScalar(field);

}

std::string MatlabInterface::getString(const std::string& name) {

  mxArray *field=getField(elem,name);
  if (!field)
    throw string("The required attribute " + name + " is missing.");

  const char *valueStr = mxArrayToString(field);
  string ret = string(valueStr);
  mxFree((void *)valueStr);
  return ret;

}

double MatlabInterface::getDouble(const std::string& name) {

  mxArray *field=getField(elem,name);
  if (!field)
    throw string("The required attribute " + name + " is missing.");
  return mxGetScalar(field);

}

double *MatlabInterface::getNativeDoubleArray(const std::string& name,std::vector<int64_t>& shape) {

  mxArray *field=getField(elem,name);
  if (!field)
    throw string("The required attribute " + name + " is missing.");

  size_t nDim = mxGetNumberOfDimensions(field);
  shape.resize(nDim);
  for(int i=0;i<nDim;i++)
    shape[i] = mxGetDimensions(field)[i];

  // Convert 1,x array to single dimension array
  if( shape[0]==1 )
    shape.erase(shape.begin());

  double *ptr = mxGetDoubles(field);
  return ptr;

}

float *MatlabInterface::getNativeFloatArray(const std::string& name,std::vector<int64_t>& shape) {

  throw string(name + ": float32 array not supported in MATLAB");

}

// --------------------------------------------------------------------------------------------------------------------

