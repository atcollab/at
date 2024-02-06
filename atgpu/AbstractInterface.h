#ifndef AT_GPU_ABSTRACTINTERFACE_H
#define AT_GPU_ABSTRACTINTERFACE_H
#include <string>
#include <vector>

// Abstract interface class
// Abstract Python, Matlab layer or CppObject debugging layer

class AbstractInterface {

public:

  // Native function
  virtual std::string getString(const std::string& name) = 0;
  virtual int getInt(const std::string& name) = 0;
  virtual double getDouble(const std::string& name) = 0;
  virtual double *getNativeDoubleArray(const std::string& name,std::vector<int64_t>& shape) = 0;
  virtual float *getNativeFloatArray(const std::string& name,std::vector<int64_t>& shape) = 0;

  // Helper function
  int getOptionalInt(const std::string& name,int defaultValue);
  double getOptionalDouble(const std::string& name,double defaultValue);
  void get1DArray(double **dest,const std::string& name,int length);
  void getOptional1DArray(double **dest,const std::string& name,int length);
  void getOptional1DArray(double **dest,const std::string& name,int *length);
  void get1DArray(float **dest,const std::string& name,int length);
  void getOptional1DArray(float **dest,const std::string& name,int length);

  // Return handle to singleton class
  static AbstractInterface *getInstance();
  static void setHandler(AbstractInterface *obj);
  static bool isValidHandler();

  // Get shape as string
  static std::string getShapeStr(std::vector<int64_t>& shape);
  static void getShapeFromStr(std::vector<int64_t>& shape,std::string& str);
  static void split(std::vector<std::string> &tokens, const std::string &text, char sep);

private:
  static AbstractInterface *handler;

};


#endif //AT_GPU_ABSTRACTINTERFACE_H
