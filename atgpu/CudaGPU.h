#ifndef AT_GPU_CUDAGPU_H
#define AT_GPU_CUDAGPU_H
#include "AbstractGPU.h"

class CudaGPU: public AbstractGPU {

public:
  CudaGPU();

  // Return list of GPU device present on the system
  std::vector<GPU_INFO> getDeviceList();

  // Return device function qualifier
  void getDeviceFunctionQualifier(std::string& ftype);

  // Compile and run the kernel
  void run(std::string& code);

  // Copy from host to device
  void hostToDevice(void *dest,void *src,size_t size);

  // Allocate device memory
  void allocDevice(void **dest,size_t size);

private:

  // Check error and throw exception in case of failure
  void cudaCheckCall(const char *funcName);
  void nvrtcCheckCall(const char *funcName);
  void addMathFunction(std::string &code);

};


#endif //AT_GPU_CUDAGPU_H
