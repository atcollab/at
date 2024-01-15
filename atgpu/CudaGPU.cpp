#include "CudaGPU.h"
#include <iostream>
#include <inttypes.h>
#include "Element.h"

using namespace std;
// -----------------------------------------------------------------------------------------------------------------

#define nvrtcCall(funcName,...) { nvrtcResult r = funcName(__VA_ARGS__);nvrtcCheckCall(#funcName,r); }
#define cudaCall(funcName,...) { CUresult r = funcName(__VA_ARGS__);cudaCheckCall(#funcName,r); }

CudaContext::CudaContext(CUDA_GPU_INFO* gpu) {

  char name[256];
  module = nullptr;
  kernel = nullptr;
  mapkernel = nullptr;
  info.name = "unknown";
  arch = "";

  cudaCall(cuDeviceGet,&cuDevice, gpu->devId);
  info = gpu->info;
  arch = gpu->arch;

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
    throw errStr;

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

  // Load the generated code and get a handle to the kernel entry point
  cudaCall(cuModuleLoadDataEx,&module, compiledCode, 0, nullptr, nullptr);
  cudaCall(cuModuleGetFunction,&kernel, module, "track");
  cudaCall(cuModuleGetFunction,&mapkernel, module, "mapbuffer");

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

// Map memory buffer
void CudaContext::mapBuffer(void **ring,uint32_t nbElement) {

  args.clear();
  args.push_back(ring);
  args.push_back(&nbElement);
  cudaCall(cuLaunchKernel, mapkernel,
           1, 1, 1,               // grid dim
           1, 1, 1,               // block dim
           0, nullptr,            // shared mem and stream
           args.data(), nullptr); // arguments

}

void CudaContext::run(uint32_t blockSize, uint64_t nbThread) {

  // Set the work item dimensions
  // Add dummy threads to allow a constant blockSize for performance
  uint32_t blockNumber = nbThread/blockSize + (((nbThread%blockSize)==0)?0:1);
  cudaCall(cuLaunchKernel, kernel,
           blockNumber, 1, 1,      // grid dim
           blockSize, 1, 1,        // block dim
           0, nullptr,             // shared mem and stream
           args.data(), nullptr);  // arguments

  // Wait end of execution
  cudaCall(cuCtxSynchronize);

}

// Check error and throw exception in case of failure
void CudaContext::cudaCheckCall(const char *funcName,CUresult r) {

  if (r != CUDA_SUCCESS) {
    const char *msg;
    cuGetErrorName(r, &msg);
    string errStr = string(funcName) + ": " + msg + " (on " + info.name + " [" + arch + "]"+ ")";
    throw errStr;
  }

}

void CudaContext::nvrtcCheckCall(const char *funcName,nvrtcResult r) {

  if(r != NVRTC_SUCCESS) {
    string errStr = string(funcName) + ": " + nvrtcGetErrorString(r) + " (on " + info.name + " [" + arch + "]"+ ")";
    throw errStr;
  }

}

// -----------------------------------------------------------------------------------------------------------------

CudaGPU::CudaGPU() {

  initErrorStr.clear();
  cudaGPUs.clear();

  // Init CUDA driver API
  CUresult status = cuInit(0);
  if( status != CUDA_SUCCESS) {
    const char *msg;
    cuGetErrorName(status, &msg);
    initErrorStr = "CUDA init failed: " + string(msg);
    return;
  }

  try {

    // Get list of CUDA devices
    char name[256];
    vector<GPU_INFO> gpuList;
    int deviceCount = 0;
    cudaCall(cuDeviceGetCount,&deviceCount);

    for(int i=0;i< deviceCount;i++) {
      CUdevice cuDevice;
      cudaCall(cuDeviceGet, &cuDevice, i);
      int min, maj, nbMP;
      cudaCall(cuDeviceGetAttribute, &min, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, cuDevice);
      cudaCall(cuDeviceGetAttribute, &maj, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, cuDevice);
      cudaCall(cuDeviceGetAttribute, &nbMP, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, cuDevice);
      cudaCall(cuDeviceGetName, name, 256, cuDevice);
      CUDA_GPU_INFO gpuInfo;
      gpuInfo.devId = i;
      gpuInfo.arch = "sm_" + to_string(maj) + to_string(min);
      gpuInfo.info.name = name;
      gpuInfo.info.version = to_string(maj) + "." + to_string(min);
      gpuInfo.info.smNumber = _ConvertSMVer2Cores(maj, min);
      gpuInfo.info.mpNumber = nbMP;
      gpuInfo.info.platform = "NVIDIA CUDA";
      cudaGPUs.push_back(gpuInfo);
    }

  } catch (string &errStr) {

    cudaGPUs.clear();
    initErrorStr = "CUDA init failed: " + errStr;

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
  if(!initErrorStr.empty())
    throw initErrorStr;
  int nbGPU = (int)cudaGPUs.size();
  if(devId>=nbGPU)
    throw string("Invalid GPU id, got " + to_string(devId) + " (Number of device: " + to_string(nbGPU) + ")");
  return new CudaContext(&cudaGPUs[devId]);
}

void CudaGPU::getDeviceFunctionQualifier(std::string& ftype) {
  ftype.assign("__device__ ");
}

void CudaGPU::getKernelFunctionQualifier(std::string& ftype) {
  // extern C to prevent from mangled name
  ftype.assign("extern \"C\" __global__ ");
}

// Return command to compute the thread id
void CudaGPU::getThreadId(std::string& command) {
  command.assign("  int threadId = blockIdx.x * blockDim.x + threadIdx.x;\n");
}

// Return global memory qualifier
void CudaGPU::getGlobalQualifier(std::string& ftype) {
  ftype.clear();
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

  if(!initErrorStr.empty())
    throw initErrorStr;

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
    info.platform = "NVIDIA CUDA";
    gpuList.push_back(info);
  }
  return gpuList;

}

// Add implementation specific function to the code
void CudaGPU::addSpecificFunctions(std::string& code) {
  if(sizeof(AT_FLOAT)==8) {
    code.append(
            "#define INF   __longlong_as_double(0x7ff0000000000000ULL)\n"
            "#define NAN   __longlong_as_double(0xfff8000000000000ULL)\n"
    );
  } else {
    code.append(
            "#define INF   __int_as_float(0x7f800000UL)\n"
            "#define NAN   __int_as_float(0x7fffffffUL)\n"
    );
  }
}

std::string CudaGPU::formatFloat(double *f) {
  char bStr[128];
  if( sizeof(AT_FLOAT)==8 ) {
#if defined(_MSC_VER) // MSVC
    sprintf(bStr, "__longlong_as_double(0x%016I64XULL)", *((uint64_t *)f));
#else
    sprintf(bStr, "__longlong_as_double(0x%" PRIx64  "ULL)", *((uint64_t *) f));
#endif
  } else {
    float f32 = (float)(*f);
    sprintf(bStr, "__int_as_float(0x%08XUL)", *((uint32_t *)&f32));
  }
  return string(bStr);
}


