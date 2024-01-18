#ifndef AT_GPU_OPENCLGPU_H
#define AT_GPU_OPENCLGPU_H
#include "AbstractGPU.h"

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 200
#define CL_HPP_TARGET_OPENCL_VERSION 200
#include <CL/opencl.hpp>

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

  // Implementation declaration
  GPU_CONTEXT_IMPL();

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

  // Implementation declaration
  GPU_IMPL();

private:

  // Check error and throw exception in case of failure
  static void openCLCheckCall(const char *funcName,cl_int r);

  std::vector<OCL_GPU_INFO> oclGPUs;

};

#endif //AT_GPU_OPENCLGPU_H
