#ifndef AT_GPU_CUDAGPU_H
#define AT_GPU_CUDAGPU_H
#include "AbstractGPU.h"
#include <cuda.h>
#include <nvrtc.h>

typedef struct {
  int devId;
  std::string arch;
  GPU_INFO info;
} CUDA_GPU_INFO;

class CudaContext final: public GPUContext {

public:

  explicit CudaContext(CUDA_GPU_INFO* gpu);
  ~CudaContext();

  // Implementation declaration
  GPU_CONTEXT_IMPL();

private:

  void cudaCheckCall(const char *funcName,CUresult r);
  void nvrtcCheckCall(const char *funcName,nvrtcResult r);

  CUcontext context;
  CUmodule module;
  CUfunction kernel;
  CUfunction mapkernel;
  std::vector<void *> args;
  GPU_INFO info;
  std::string arch;
  CUdevice cuDevice;

};

class CudaGPU: public AbstractGPU {

public:

  CudaGPU();

  // Implementation declaration
  GPU_IMPL();

private:

  // Check error and throw exception in case of failure
  static void cudaCheckCall(const char *funcName,CUresult r);

  std::vector<CUDA_GPU_INFO> cudaGPUs;

};


#endif //AT_GPU_CUDAGPU_H
