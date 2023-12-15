#ifndef AT_GPU_CUDAGPU_H
#define AT_GPU_CUDAGPU_H
#include "AbstractGPU.h"
#include <cuda.h>
#include <nvrtc.h>

class CudaContext final: public GPUContext {

public:

  explicit CudaContext(int devId) noexcept;
  ~CudaContext() final;

  // Empty kernel parameter
  void resetArg() final;

  // Set kernel parameter
  void addArg(size_t argSize,void *value) final;

  // Run the kernel
  void run(uint32_t gridSize,uint64_t nbThread) final;

  // Compile and load the kernel
  void compile(std::string& code) final;

  // Copy from host to device
  void hostToDevice(void *dest,void *src,size_t size) final;

  // Copy from device to host
  void deviceToHost(void *dest,void *src,size_t size) final;

  // Allocate device memory
  void allocDevice(void **dest,size_t size) final;

  // Free device memory
  void freeDevice(void *dest) final;

  // Return device name
  std::string& name();

private:

  void cudaCheckCall(const char *funcName,CUresult r);
  void nvrtcCheckCall(const char *funcName,nvrtcResult r);

  CUcontext context;
  CUmodule module;
  CUfunction kernel;
  int devId;
  std::vector<void *> args;
  GPU_INFO info;
  std::string arch;
  CUdevice cuDevice;

};

class CudaGPU: public AbstractGPU {

public:
  CudaGPU();

  // Return list of GPU device present on the system
  std::vector<GPU_INFO> getDeviceList() override;

  // Create a context for the given GPU
  GPUContext *createContext(int devId) override;

  // Return device function qualifier
  void getDeviceFunctionQualifier(std::string& ftype) override;

  // Check error and throw exception in case of failure
  static void cudaCheckCall(const char *funcName,CUresult r);
  static int _ConvertSMVer2Cores(int major,int minor);

};


#endif //AT_GPU_CUDAGPU_H
