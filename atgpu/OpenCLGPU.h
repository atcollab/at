#ifndef AT_GPU_OPENCLGPU_H
#define AT_GPU_OPENCLGPU_H
#include "AbstractGPU.h"

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 200
#define CL_HPP_TARGET_OPENCL_VERSION 200
#if defined(_MSC_VER) // MSVC
#include <CL/opencl.h>
#else
#include <CL/opencl.hpp>
#endif

typedef struct {
  cl_platform_id platform_id;
  cl_device_id device_id;
  bool fp64;
  GPU_INFO info;
} OCL_GPU_INFO;

typedef struct {
  size_t size;
  void *arg;
} ARG;

class OpenCLContext final: public GPUContext {

public:

  explicit OpenCLContext(OCL_GPU_INFO *gpu);
  ~OpenCLContext() final;

  // Empty kernel parameter
  void resetArg() final;

  // Set kernel parameter
  void addArg(size_t argSize,void *value) final;

  // Run the kernel
  void run(uint32_t blockSize, uint64_t nbThread) final;

  // Compile and load the kernel
  void compile(std::string& code) final;

  // Map memory buffer
  void mapBuffer(void **ring,uint32_t nbElement);

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

  // Return number fo core
  int coreNumber();

private:

  void openCLCheckCall(const char *funcName,cl_int r);

  std::vector<ARG> args;
  GPU_INFO info;
  cl_platform_id platform_id;
  cl_device_id device_id;
  cl_context context;
  cl_command_queue commands;
  cl_program program;
  cl_kernel kernel;
  cl_kernel mapkernel;

};

class OpenCLGPU: public AbstractGPU {

public:
  OpenCLGPU();

  // Return list of GPU device present on the system
  std::vector<GPU_INFO> getDeviceList() override;

  // Create a context for the given GPU
  GPUContext *createContext(int devId) override;

  // Return device function qualifier
  void getDeviceFunctionQualifier(std::string& ftype) override;

  // Return kernel function qualifier
  void getKernelFunctionQualifier(std::string& ftype);

  // Return global memory qualifier
  void getGlobalQualifier(std::string& ftype);

  // Add implementation specific function to the code
  void addSpecificFunctions(std::string& code);

  // Return command to compute the thread id
  void getThreadId(std::string& command);

  // Format a double
  std::string formatFloat(double *f);

private:

  // Check error and throw exception in case of failure
  static void openCLCheckCall(const char *funcName,cl_int r);

  std::vector<OCL_GPU_INFO> oclGPUs;

};

#endif //AT_GPU_OPENCLGPU_H
