#include "CudaGPU.h"
#include <iostream>

using namespace std;
// -----------------------------------------------------------------------------------------------------------------

#define nvrtcCall(funcName,...) { nvrtcResult r = funcName(__VA_ARGS__);nvrtcCheckCall(#funcName,r); }
#define cudaCall(funcName,...) { CUresult r = funcName(__VA_ARGS__);cudaCheckCall(#funcName,r); }

CudaContext::CudaContext(int devId) noexcept: devId(devId)  {

  char name[256];
  module = nullptr;
  kernel = nullptr;

  cudaCall(cuDeviceGet,&cuDevice, devId);
  int min,maj,nbMP;
  cudaCall(cuDeviceGetAttribute,&min,CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR,cuDevice);
  cudaCall(cuDeviceGetAttribute,&maj,CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR,cuDevice);
  cudaCall(cuDeviceGetAttribute,&nbMP,CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT,cuDevice);
  cudaCall(cuDeviceGetName, name, 256, cuDevice );
  info.name = name;
  info.version = to_string(maj) + "." + to_string(min);
  info.mpNumber = nbMP;
  info.smNumber = CudaGPU::_ConvertSMVer2Cores(maj,min);
  arch = "sm_" + to_string(maj) + to_string(min);

  //cudaCall(cuDevicePrimaryCtxRetain, &context, cuDevice);
  cudaCall(cuCtxCreate,&context, CU_CTX_SCHED_SPIN, cuDevice);
  cudaCall(cuCtxSetSharedMemConfig, CU_SHARED_MEM_CONFIG_EIGHT_BYTE_BANK_SIZE);
  cudaCall(cuCtxSetCacheConfig , CU_FUNC_CACHE_PREFER_L1);

}

int CudaContext::coreNumber() {
  return info.mpNumber * info.smNumber;
}

CudaContext::~CudaContext() {
  if(module) cuModuleUnload(module);
  cuCtxDestroy(context);
  //cudaCall(cuDevicePrimaryCtxRelease, cuDevice);
}

void CudaContext::hostToDevice(void *dest,void *src,size_t size) {
  cudaCall(cuMemcpyHtoD,(CUdeviceptr)dest,src,size);
}

void CudaContext::deviceToHost(void *dest,void *src,size_t size) {
  cudaCall(cuMemcpyDtoH,dest,(CUdeviceptr)src,size);
}

void CudaContext::allocDevice(void **dest,size_t size) {
  cudaCall(cuMemAlloc,(CUdeviceptr *)dest,size);
}

void CudaContext::freeDevice(void *dest) {
  cudaCall(cuMemFree,(CUdeviceptr)dest);
}

void CudaContext::compile(string& code) {

  // Create the program
  nvrtcProgram prog;
  nvrtcCall(nvrtcCreateProgram,
            &prog,         // prog
            code.c_str(),  // buffer
            "track.cu",    // name
            0,             // numHeaders
            nullptr,       // headers
            nullptr);      // includeNames

  // Compile the program
  string archOpt = "-arch=" + arch;
  const char *opts[] = {archOpt.c_str()};
  nvrtcResult compileResult = nvrtcCompileProgram(prog,  // prog
                                                  1,     // numOptions
                                                  opts); // options

  if (compileResult != NVRTC_SUCCESS) {

    AbstractGPU::outputCode(code);

    size_t logSize;
    nvrtcCall(nvrtcGetProgramLogSize,prog, &logSize);
    char *log = new char[logSize];
    nvrtcCall(nvrtcGetProgramLog,prog, log);
    std::cout << log << '\n';
    delete[] log;

    // Obtain compilation log from the program.

    nvrtcDestroyProgram(&prog);
    string errStr = "nvrtcCompileProgram failed: " + string(nvrtcGetErrorString(compileResult));
    throw(errStr);

  }

#if 0
  // nvcc --cubin --ptxas-options="-m64 -arch=sm_86 --verbose -O3" code.ptx
  char *compiledCode = new char[36656];
  FILE *f = fopen("code.cubin","r");
  fread(compiledCode,1,36656,f);
  fclose(f);
#endif

#if 0
  // Obtain PTX from the program.
  size_t ptxSize;
  nvrtcCall(nvrtcGetPTXSize,prog, &ptxSize);
  char *compiledCode = new char[ptxSize];
  nvrtcCall(nvrtcGetPTX,prog,compiledCode);
  FILE *f = fopen("code.ptx","w");
  fwrite(compiledCode,1,ptxSize,f);
  fclose(f);
#endif

#if 1
  // Obtain BIN from the program.
  size_t binSize;
  nvrtcCall(nvrtcGetCUBINSize,prog, &binSize);
  char *compiledCode = new char[binSize];
  nvrtcCall(nvrtcGetCUBIN,prog, compiledCode);
#endif

  // Destroy the program.
  nvrtcCall(nvrtcDestroyProgram,&prog);

  // Load the generated code and get a handle to the kernel.
  cudaCall(cuModuleLoadDataEx,&module, compiledCode, 0, nullptr, nullptr);

}

std::string& CudaContext::name() {
  return info.name;
}

void CudaContext::resetArg() {
  args.clear();
}

void CudaContext::addArg(size_t argSize,void *value) {
  args.push_back(value);
}

void CudaContext::run(uint32_t blockSize, uint64_t nbThread) {

  cudaCall(cuModuleGetFunction,&kernel, module, "track");

  if(nbThread < blockSize ) {

    cudaCall(cuLaunchKernel, kernel,
             1, 1, 1,               // grid dim
             blockSize, 1, 1,       // block dim
             0, nullptr,            // shared mem and stream
             args.data(), nullptr); // arguments

  } else {

    if(nbThread % blockSize != 0  ) {
      // TODO: Handle this
      throw string("nbThread (" + to_string(nbThread) + ") must be a multiple of GPU_BLOCK_SIZE (" + to_string(blockSize) + ")");
    }

    cudaCall(cuLaunchKernel, kernel,
             nbThread / blockSize, 1, 1,  // grid dim
             blockSize, 1, 1,           // block dim
             0, nullptr,               // shared mem and stream
             args.data(), nullptr);    // arguments

  }

  // Wait end of execution
  cudaCall(cuCtxSynchronize);

}

// Check error and throw exception in case of failure
void CudaContext::cudaCheckCall(const char *funcName,CUresult r) {

  if (r != CUDA_SUCCESS) {
    const char *msg;
    cuGetErrorName(r, &msg);
    string errStr = info.name + "[" + info.version + "] " + string(funcName) + ": " + msg;
    throw errStr;
  }

}

void CudaContext::nvrtcCheckCall(const char *funcName,nvrtcResult r) {

  if(r != NVRTC_SUCCESS) {
    string errStr = info.name + "[" + info.version + "] " + string(funcName) + ": " + nvrtcGetErrorString(r);
    throw errStr;
  }

}

// -----------------------------------------------------------------------------------------------------------------

CudaGPU::CudaGPU() {

  CUresult r = cuInit(0);
  if( r != CUDA_SUCCESS) {
    const char *msg;
    cuGetErrorName(r, &msg);
    string errStr = "Cannot initialise CUDA:" + string(msg);
    throw errStr;
  }

};

void CudaGPU::cudaCheckCall(const char *funcName,CUresult r) {

  if (r != CUDA_SUCCESS) {
    const char *msg;
    cuGetErrorName(r, &msg);
    string errStr = string(funcName) + ": " + msg;
    throw errStr;
  }

}

GPUContext *CudaGPU::createContext(int devId) {
  return new CudaContext(devId);
}

void CudaGPU::getDeviceFunctionQualifier(std::string& ftype) {
  ftype.assign("__device__");
}

// Return number of stream processors
int CudaGPU::_ConvertSMVer2Cores(int major,int minor) {

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


// Return list of GPU device present on the system
std::vector<GPU_INFO> CudaGPU::getDeviceList() {

  char name[256];
  vector<GPU_INFO> gpuList;
  int deviceCount = 0;
  cudaCall(cuDeviceGetCount,&deviceCount);

  for(int i=0;i< deviceCount;i++) {
    CUdevice cuDevice;
    cudaCall(cuDeviceGet,&cuDevice, i);
    int min,maj,nbMP;
    cudaCall(cuDeviceGetAttribute,&min,CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR,cuDevice);
    cudaCall(cuDeviceGetAttribute,&maj,CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR,cuDevice);
    cudaCall(cuDeviceGetAttribute,&nbMP,CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT,cuDevice);
    cudaCall(cuDeviceGetName, name, 256, cuDevice );
    GPU_INFO info;
    info.name = name;
    info.version = to_string(maj)+ "." + to_string(min);
    info.smNumber = _ConvertSMVer2Cores(maj,min);
    info.mpNumber = nbMP;
    gpuList.push_back(info);
  }
  return gpuList;

}


