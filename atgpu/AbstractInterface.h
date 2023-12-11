#ifndef AT_GPU_ABSTRACTINTERFACE_H
#define AT_GPU_ABSTRACTINTERFACE_H
#include <string>
#include <vector>
#include "Element.h"

// Abstract interface class
// Abstract Python or Matlab layer
class AbstractInterface {

public:

  // Native function
  virtual std::string getString(const std::string& name) = 0;
  virtual int getInt(const std::string& name) = 0;
  virtual AT_FLOAT getDouble(const std::string& name) = 0;
  virtual AT_FLOAT *getNativeDoubleArray(const std::string& name,std::vector<int64_t>& shape) = 0;

  // Helper function
  int getOptionalInt(const std::string& name,int defaultValue);
  AT_FLOAT getOptionalDouble(const std::string& name,AT_FLOAT defaultValue);
  AT_FLOAT *getDoubleArray(const std::string& name,std::vector<int64_t> expectedShape);
  AT_FLOAT *getOptionalDoubleArray(const std::string& name,std::vector<int64_t> expectedShape);
  void get1DArray(AT_FLOAT **dest,const std::string& name,int length);
  void getOptional1DArray(AT_FLOAT **dest,const std::string& name,int length);

  // Return handle to singleton class
  static AbstractInterface *getInstance();
  static void setHandler(AbstractInterface *obj);

  // Get shape as string
  static std::string getShapeStr(std::vector<int64_t>& shape);

private:
  static AbstractInterface *handler;

};


#endif //AT_GPU_ABSTRACTINTERFACE_H
