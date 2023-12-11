#include "CudaGPU.h"
#include <cuda_runtime.h>
#include <nvrtc.h>
#include <iostream>

using namespace std;

#define cudaCall(funcName,...) funcName(__VA_ARGS__);cudaCheckCall(#funcName);
#define nvrtcCall(funcName,...) funcName(__VA_ARGS__);nvrtcCheckCall(#funcName);

CudaGPU::CudaGPU() {

}

// Copy from host to device
void CudaGPU::hostToDevice(void *dest,void *src,size_t size) {
  cudaCall(cudaMemcpy,dest,src,size,cudaMemcpyHostToDevice);
}

void CudaGPU::getDeviceFunctionQualifier(std::string& ftype) {
  ftype.assign("__device__");
}

// Allocate device memory
void CudaGPU::allocDevice(void **dest,size_t size) {
  cudaCall(cudaMalloc,dest,size);
}

// Return number of stream processors
static int _ConvertSMVer2Cores(int major,int minor) {

  // Defines for GPU Architecture types (using the SM version to determine
  // the # of cores per SM
  typedef struct {
    int SM;  // 0xMm (hexidecimal notation), M = SM Major version,
    // and m = SM minor version
    int Cores;
  } sSMtoCores;

  sSMtoCores nGpuArchCoresPerSM[] = {
          {0x20, 32},  // Fermi
          {0x21, 48},  // Fermi
          {0x30, 192}, // Kepler
          {0x32, 192}, // Kepler
          {0x35, 192}, // Kepler
          {0x37, 192}, // Kepler
          {0x50, 128}, // Maxwell
          {0x52, 128}, // Maxwell
          {0x53, 128}, // Maxwell
          {0x60,  64}, // Pascal
          {0x61, 128}, // Volta
          {0x62, 128}, // Volta
          {0x70,  64}, // Turing
          {0x72,  64},
          {0x75,  64}, // Turing
          {0x80,  64}, // Ampere
          {0x86, 128}, // Ampere
          {0x87, 128}, // Ampere
          {0x89, 128}, // Ampere
          {0x90, 128}, // Hopper
          {-1, -1} };

  int index = 0;

  while(nGpuArchCoresPerSM[index].SM != -1) {
    if(nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
      return nGpuArchCoresPerSM[index].Cores;
    }

    index++;
  }

  return 0;

}

// Check error and throw exception in case of failure
void CudaGPU::cudaCheckCall(const char *funcName) {

  cudaError_t errCode=cudaGetLastError();
  if(errCode != cudaSuccess) {
    string errStr = string(funcName) + ": " + cudaGetErrorString(errCode);
    throw errStr;
  }

}

void CudaGPU::nvrtcCheckCall(const char *funcName) {

  cudaError_t errCode=cudaGetLastError();
  if(errCode != cudaSuccess) {
    string errStr = string(funcName) + ": " + cudaGetErrorString(errCode);
    throw errStr;
  }

}

// Return list of GPU device present on the system
std::vector<GPU_INFO> CudaGPU::getDeviceList() {

  vector<GPU_INFO> gpuList;
  int deviceCount = 0;
  cudaCall(cudaGetDeviceCount,&deviceCount);

  for(int i=0;i< deviceCount;i++) {

    cudaDeviceProp deviceProp;
    cudaCall(cudaGetDeviceProperties,&deviceProp,i);
    GPU_INFO info;
    info.name = deviceProp.name;
    info.version = to_string(deviceProp.major)+ "." + to_string(deviceProp.minor);
    info.smNumber = _ConvertSMVer2Cores(deviceProp.major,deviceProp.minor);
    info.mpNumber = deviceProp.multiProcessorCount;
    gpuList.push_back(info);
  }
  return gpuList;

}

void CudaGPU::run(string& code) {

  // Create the program
  nvrtcProgram prog;
  nvrtcCall(nvrtcCreateProgram,
                     &prog,         // prog
                     code.c_str(),  // buffer
                     "track.cu",    // name
                     0,             // numHeaders
                     nullptr,       // headers
                     nullptr);      // includeNames

  // Compile the program with fmad disabled.
  // Note: Can specify GPU target architecture explicitly with '-arch' flag.
  const char *opts[] = {"--gpu-architecture=compute_61"};
  nvrtcResult compileResult = nvrtcCompileProgram(prog,  // prog
                                                  1,     // numOptions
                                                  opts); // options
  if (compileResult != NVRTC_SUCCESS) {

    outputCode(code);

    // Obtain compilation log from the program.
    size_t logSize;
    nvrtcCall(nvrtcGetProgramLogSize,prog, &logSize);
    char *log = new char[logSize];
    nvrtcCall(nvrtcGetProgramLog,prog, log);
    std::cout << log << '\n';
    delete[] log;

    nvrtcDestroyProgram(&prog);
    string errStr = "nvrtcCompileProgram failed: " + string(nvrtcGetErrorString(compileResult));
    throw(errStr);

  }

  /*

  // Obtain PTX from the program.
  size_t ptxSize;
  NVRTC_SAFE_CALL(nvrtcGetPTXSize(prog, &ptxSize));
  char *ptx = new char[ptxSize];
  NVRTC_SAFE_CALL(nvrtcGetPTX(prog, ptx));
  // Destroy the program.
  NVRTC_SAFE_CALL(nvrtcDestroyProgram(&prog));
  // Load the generated PTX and get a handle to the SAXPY kernel.
  CUdevice cuDevice;
  CUcontext context;
  CUmodule module;
  CUfunction kernel;
  CUDA_SAFE_CALL(cuInit(0));
  CUDA_SAFE_CALL(cuDeviceGet(&cuDevice, 0));
  CUDA_SAFE_CALL(cuCtxCreate(&context, 0, cuDevice));
  CUDA_SAFE_CALL(cuModuleLoadDataEx(&module, ptx, 0, 0, 0));
  CUDA_SAFE_CALL(cuModuleGetFunction(&kernel, module, "track"));
  */

}

