#include "OpenCLGPU.h"
#include <iostream>

// Maximum number of GPU on a host
#define MAX_OCL_GPU 16
// Maximum number of platform on a host
#define MAX_OCL_PLATFORM 16

using namespace std;

#ifdef __APPLE__
#define clCreateCommandQueueWithProperties(ctxt,devid,xx,err) clCreateCommandQueue(ctxt,devid,0,err)
#define OPTIONS ""
#else
#define OPTIONS "-D__GNUC__"
#endif /*__APPLE__*/

#define openCLCall(funcName,...) { cl_int r = funcName(__VA_ARGS__);openCLCheckCall(#funcName,r); }

static string getCLErrorString(cl_int r) {
  string msg;
  switch (r) {
    case CL_DEVICE_NOT_FOUND:                   msg = "Device not found.";break;
    case CL_DEVICE_NOT_AVAILABLE:               msg = "Device not available";break;
    case CL_COMPILER_NOT_AVAILABLE:             msg = "Compiler not available";break;
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:      msg = "Memory object allocation failure";break;
    case CL_OUT_OF_RESOURCES:                   msg = "Out of resources";break;
    case CL_OUT_OF_HOST_MEMORY:                 msg = "Out of host memory";break;
    case CL_PROFILING_INFO_NOT_AVAILABLE:       msg = "Profiling information not available";break;
    case CL_MEM_COPY_OVERLAP:                   msg = "Memory copy overlap";break;
    case CL_IMAGE_FORMAT_MISMATCH:              msg = "Image format mismatch";break;
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:         msg = "Image format not supported";break;
    case CL_BUILD_PROGRAM_FAILURE:              msg = "Program build failure";break;
    case CL_MAP_FAILURE:                        msg = "Map failure";break;
    case CL_INVALID_VALUE:                      msg = "Invalid value";break;
    case CL_INVALID_DEVICE_TYPE:                msg = "Invalid device type";break;
    case CL_INVALID_PLATFORM:                   msg = "Invalid platform";break;
    case CL_INVALID_DEVICE:                     msg = "Invalid device";break;
    case CL_INVALID_CONTEXT:                    msg = "Invalid context";break;
    case CL_INVALID_QUEUE_PROPERTIES:           msg = "Invalid queue properties";break;
    case CL_INVALID_COMMAND_QUEUE:              msg = "Invalid command queue";break;
    case CL_INVALID_HOST_PTR:                   msg = "Invalid host pointer";break;
    case CL_INVALID_MEM_OBJECT:                 msg = "Invalid memory object";break;
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    msg = "Invalid image format descriptor";break;
    case CL_INVALID_IMAGE_SIZE:                 msg = "Invalid image size";break;
    case CL_INVALID_SAMPLER:                    msg = "Invalid sampler";break;
    case CL_INVALID_BINARY:                     msg = "Invalid binary";break;
    case CL_INVALID_BUILD_OPTIONS:              msg = "Invalid build options";break;
    case CL_INVALID_PROGRAM:                    msg = "Invalid program";break;
    case CL_INVALID_PROGRAM_EXECUTABLE:         msg = "Invalid program executable";break;
    case CL_INVALID_KERNEL_NAME:                msg = "Invalid kernel name";break;
    case CL_INVALID_KERNEL_DEFINITION:          msg = "Invalid kernel definition";break;
    case CL_INVALID_KERNEL:                     msg = "Invalid kernel";break;
    case CL_INVALID_ARG_INDEX:                  msg = "Invalid argument index";break;
    case CL_INVALID_ARG_VALUE:                  msg = "Invalid argument value";break;
    case CL_INVALID_ARG_SIZE:                   msg = "Invalid argument size";break;
    case CL_INVALID_KERNEL_ARGS:                msg = "Invalid kernel arguments";break;
    case CL_INVALID_WORK_DIMENSION:             msg = "Invalid work dimension";break;
    case CL_INVALID_WORK_GROUP_SIZE:            msg = "Invalid work group size";break;
    case CL_INVALID_WORK_ITEM_SIZE:             msg = "Invalid work item size";break;
    case CL_INVALID_GLOBAL_OFFSET:              msg = "Invalid global offset";break;
    case CL_INVALID_EVENT_WAIT_LIST:            msg = "Invalid event wait list";break;
    case CL_INVALID_EVENT:                      msg = "Invalid event";break;
    case CL_INVALID_OPERATION:                  msg = "Invalid operation";break;
    case CL_INVALID_GL_OBJECT:                  msg = "Invalid OpenGL object";break;
    case CL_INVALID_BUFFER_SIZE:                msg = "Invalid buffer size";break;
    case CL_INVALID_MIP_LEVEL:                  msg = "Invalid mip-map level";break;
    default: msg = "Unknown error code #" + to_string(r);
  }
  return msg;
}

// -----------------------------------------------------------------------------------------------------------------

OpenCLContext::OpenCLContext(OCL_GPU_INFO *gpu)  {

  info = gpu->info;

  cl_int err;
  platform_id = gpu->platform_id;
  device_id = gpu->device_id;
  cl_context_properties cps[] = {CL_CONTEXT_PLATFORM,(cl_context_properties)platform_id,0};
  context = clCreateContext(cps, 1, &device_id, nullptr, nullptr, &err);
  if(!context)
    openCLCheckCall("clCreateContext",err);

  commands = clCreateCommandQueueWithProperties(context, device_id, nullptr, &err);
  if(!commands)
    openCLCheckCall("clCreateCommandQueue",err);

  program = nullptr;
  kernel = nullptr;
  mapkernel = nullptr;

}

int OpenCLContext::coreNumber() {
  return info.mpNumber;
}

OpenCLContext::~OpenCLContext() {
  if(mapkernel) clReleaseKernel(mapkernel);
  if(kernel) clReleaseKernel(kernel);
  if(program) clReleaseProgram(program);
  clReleaseCommandQueue(commands);
  clReleaseContext(context);
}

void OpenCLContext::hostToDevice(void *dest,void *src,size_t size) {
  openCLCall(clEnqueueWriteBuffer,commands, (cl_mem)dest, CL_TRUE, 0, size, src, 0, nullptr, nullptr);
}

void OpenCLContext::deviceToHost(void *dest,void *src,size_t size) {
  openCLCall(clEnqueueReadBuffer,commands, (cl_mem)src, CL_TRUE, 0, size, dest, 0, nullptr, nullptr);
}

void OpenCLContext::allocDevice(void **dest,size_t size) {
  cl_int err;
  *dest = clCreateBuffer(context, CL_MEM_READ_WRITE,size,nullptr,&err);
  openCLCheckCall("clCreateBuffer",err);
}

void OpenCLContext::freeDevice(void *dest) {
  openCLCall(clReleaseMemObject,(cl_mem)dest);
}

void OpenCLContext::compile(string& code) {

  cl_int err;
  const char *progSrc = code.c_str();
  program = clCreateProgramWithSource(context, 1, &progSrc, nullptr, &err);
  openCLCheckCall("clCreateCommandQueue",err);
  // const char *opts = "-D__GNUC__"; // -D__GNUC__ for structure alignment (same directive as gcc)
  const char *opts = OPTIONS;
  err = clBuildProgram(program, 1, &device_id, opts, nullptr, nullptr);
  if(err < 0) {

    AbstractGPU::outputCode(code);

    size_t logSize = 0;
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, nullptr, &logSize);
    if( logSize>0 ) {
      char *log = new char[logSize + 1];
      clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, logSize + 1, log, nullptr);
      std::cout << log << std::endl;
      delete[] log;
    }

    clReleaseProgram(program);
    program = nullptr;
    string errStr = "clBuildProgram failed: " + getCLErrorString(err);
    throw errStr;

  }

  kernel = clCreateKernel(program, "track", &err);
  if(err < 0) {
    clReleaseProgram(program);
    program = nullptr;
    string errStr = "clCreateKernel failed: " + getCLErrorString(err);
    throw errStr;
  }

  mapkernel = clCreateKernel(program, "mapbuffer", &err);
  if(err < 0) {
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    program = nullptr;
    kernel = nullptr;
    string errStr = "clCreateKernel failed: " + getCLErrorString(err);
    throw errStr;
  }

}

std::string& OpenCLContext::name() {
  return info.name;
}

void OpenCLContext::resetArg() {
  args.clear();
}

void OpenCLContext::addArg(size_t argSize,void *value) {
  ARG a;
  a.size = argSize;
  a.arg = value;
  args.push_back(a);
}

// Map memory buffer
void OpenCLContext::mapBuffer(void **ring,uint32_t nbElement) {

  openCLCall(clSetKernelArg,mapkernel,0,sizeof(void *), ring);
  openCLCall(clSetKernelArg,mapkernel,1,sizeof(uint32_t), &nbElement);
  size_t globalSize[] = {1};
  openCLCall(clEnqueueNDRangeKernel,commands, mapkernel, 1, nullptr, globalSize, nullptr, 0, nullptr, nullptr);
  clFinish(commands);

}

void OpenCLContext::run(uint64_t nbThread) {

  //set the kernel arguments
  cl_int err;
  for(size_t i=0;i<args.size();i++) {
    err = clSetKernelArg( kernel, (cl_uint)i, (cl_uint)args[i].size, args[i].arg);
    if( err<0 ) {
      string errStr = "Argument #" + to_string(i) + " size:" + to_string(args[i].size) + " " + getCLErrorString(err);
      throw errStr;
    }
  }

  uint32_t blockSize;
  uint32_t blockNumber;

  if( nbThread<coreNumber() ) {
    // Launch 1 thread per core (here, a core is a multiprocessor)
    blockSize = 1;
    blockNumber = nbThread;
  } else {
    // Add dummy threads to allow a constant blockSize for performance
    // Choose 64 (2 warps) which seems a good compromise
    blockSize = 64;
    blockNumber = nbThread / blockSize + (((nbThread % blockSize) == 0) ? 0 : 1);
  }

  size_t globalSize[] = {blockSize, blockNumber};
  openCLCall(clEnqueueNDRangeKernel, commands, kernel, 2, nullptr, globalSize, nullptr, 0, nullptr, nullptr);
  clFinish(commands);

}

// Check error and throw exception in case of failure
void OpenCLContext::openCLCheckCall(const char *funcName,int r) {

  if (r != CL_SUCCESS) {
    string errStr = string(funcName) + ": " + getCLErrorString(r) + " (on " + info.name + ")";
    throw errStr;
  }

}

// -----------------------------------------------------------------------------------------------------------------

OpenCLGPU::OpenCLGPU() {

  // Detect all OpenCL capable devices
  oclGPUs.clear();
  initErrorStr.clear();
  implementationStr = "OpenCL";

  try {

    cl_uint plalformCount;
    cl_platform_id platform[MAX_OCL_PLATFORM];
    openCLCall(clGetPlatformIDs, MAX_OCL_PLATFORM, platform, &plalformCount);
    vector<GPU_INFO> gpuList;

    for (cl_uint p = 0; p < plalformCount; p++) {

      char pname[256];
      char pversion[256];
      openCLCall(clGetPlatformInfo, platform[p], CL_PLATFORM_NAME, sizeof(pname), pname, nullptr);
      openCLCall(clGetPlatformInfo, platform[p], CL_PLATFORM_VERSION, sizeof(pversion), pversion, nullptr);

      try {

        cl_device_id device_id[MAX_OCL_GPU];
        cl_uint deviceCount;
        openCLCall(clGetDeviceIDs, platform[p], CL_DEVICE_TYPE_CPU | CL_DEVICE_TYPE_GPU, MAX_OCL_GPU, device_id, &deviceCount);

        for (cl_uint i = 0; i < deviceCount; i++) {

          char name[256];
          cl_uint mpNumer;
          cl_uint maj=0;
          cl_uint min=0;
          char clExtensions[16384];

          openCLCall(clGetDeviceInfo, device_id[i], CL_DEVICE_NAME, sizeof(name), name, nullptr);
          openCLCall(clGetDeviceInfo, device_id[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(mpNumer), &mpNumer, nullptr);
          openCLCall(clGetDeviceInfo, device_id[i], CL_DEVICE_EXTENSIONS, sizeof(clExtensions), &clExtensions, nullptr);

#ifdef CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
          try {
            // Retrieve CUDA capabilities
            openCLCall(clGetDeviceInfo, device_id[i], CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof(maj), &maj, nullptr);
            openCLCall(clGetDeviceInfo, device_id[i], CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV, sizeof(min), &min, nullptr);
          } catch (string& errStr) {
            // Likely not a CUDA device
          }
#endif

          OCL_GPU_INFO oInfo;
          oInfo.device_id = device_id[i];
          oInfo.platform_id = platform[p];
          oInfo.info.name = name;
          oInfo.info.version = to_string(maj) + "." + to_string(min);
          oInfo.info.mpNumber = mpNumer;
          oInfo.info.platform = string(pname) + " " + string(pversion);
          string clExtStr = clExtensions;
          oInfo.fp64 = clExtStr.find("_fp64") != string::npos;
          oclGPUs.push_back(oInfo);

        }

      } catch (string& errStr) {
        cerr << "Warning, Platform " << pname << ": " << errStr << endl;
      }

    }
  } catch (string &errStr) {

    oclGPUs.clear();
    initErrorStr = "OpenCL init failed: " + errStr;

  }

};

void OpenCLGPU::openCLCheckCall(const char *funcName,int r) {

  if (r != CL_SUCCESS) {
    string errStr = string(funcName) + ": " + getCLErrorString(r);
    throw errStr;
  }

}

GPUContext *OpenCLGPU::createContext(int devId) {
  if(!initErrorStr.empty())
    throw initErrorStr;
  int nbGPU = (int)oclGPUs.size();
  if(devId>=nbGPU)
    throw string("Invalid GPU id, got " + to_string(devId) + " (Number of device: " + to_string(nbGPU) + ")");
  if( !oclGPUs[devId].fp64 )
    throw string(oclGPUs[devId].info.name + " has no FP64 support");
  return new OpenCLContext(&oclGPUs[devId]);
}

void OpenCLGPU::getDeviceFunctionQualifier(std::string& ftype) {
  ftype.clear();
}

void OpenCLGPU::getKernelFunctionQualifier(std::string& ftype) {
  ftype.assign("__kernel ");
}

// Return global memory qualifier
void OpenCLGPU::getGlobalQualifier(std::string& ftype) {
  ftype.assign("__global ");
}

// Return command to compute the thread id
void OpenCLGPU::getThreadId(std::string& command) {
  command.assign("  int threadId = get_global_id(0) + get_global_size(0) * get_global_id(1);\n");
}

// Return list of GPU device present on the system
vector<GPU_INFO> OpenCLGPU::getDeviceList() {
  if(!initErrorStr.empty())
    throw initErrorStr;
  vector<GPU_INFO> gpuList;
  for(auto & oclGPU : oclGPUs)
    gpuList.push_back(oclGPU.info);
  return gpuList;
}

// Add implementation specific function to the code
void OpenCLGPU::addSpecificFunctions(std::string& code) {
  code.append("#define INF INFINITY\n");
}

std::string OpenCLGPU::formatFloat(double *f) {
  char bStr[128];
  sprintf(bStr, "(AT_FLOAT)%.16f", *f);
  return string(bStr);
}

