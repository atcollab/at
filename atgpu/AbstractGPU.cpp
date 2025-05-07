// Abstract class for GPU context
// Abstract OpenCL or CUDA API

#if defined(CUDA)
#include "CudaGPU.h"
#elif defined(OPENCL)
#include "OpenCLGPU.h"
#else
#error "Please provide a definition for GPU implementation (CUDA or OPENCL) !"
#endif

#include "iostream"

// Timing definition for profiling
#if defined(_MSC_VER)

#include <windows.h>
LARGE_INTEGER perfTickStart;
double perfTicksPerSec;
LARGE_INTEGER qwTicksPerSec;

#else

#include <sys/time.h>
time_t tickStart;

#endif

#include "AbstractInterface.h"

using namespace std;

// Create GPU implementation object
AbstractGPU *AbstractGPU::gpuHandler = nullptr;
AbstractGPU *AbstractGPU::getInstance() {
  if( gpuHandler== nullptr ) {
    initTimer();
#if defined(CUDA)
    gpuHandler = new CudaGPU();
#elif defined(OPENCL)
    gpuHandler = new OpenCLGPU();
#endif
  }
  return gpuHandler;
}

// Initialise timer
void AbstractGPU::initTimer() {

#if defined(_MSC_VER)
  QueryPerformanceFrequency(&qwTicksPerSec);
  QueryPerformanceCounter(&perfTickStart);
  perfTicksPerSec = (double)qwTicksPerSec.QuadPart;
#else
  tickStart = time(NULL);
#endif

}

// Return number of seconds
double AbstractGPU::get_ticks() {

#if defined(_MSC_VER)
  LARGE_INTEGER t,dt;
  QueryPerformanceCounter(&t);
  dt.QuadPart = t.QuadPart - perfTickStart.QuadPart;
  return (double)(dt.QuadPart) / perfTicksPerSec;
#else
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (double)(tv.tv_sec - tickStart) + (double)tv.tv_usec / 1e6;
#endif

}

// Output code with line number
void AbstractGPU::outputCode(std::string& code) {

  vector<string> lines;
  AbstractInterface::split(lines,code,'\n');
  for(size_t i=0;i<lines.size();i++) {
    char tmp[256];
    sprintf(tmp,"%04d: ",(int)(i+1));
    cout << tmp << lines[i] << endl;
  }

}

// Return implementation name
std::string& AbstractGPU::implName() {
  return implementationStr;
}

void AbstractGPU::addUtilsFunctions(std::string &code) {

  // Macro used for mapping buffer in GPU memory space (see Lattice::mapBuffers)
  code.append("#define MAP_FLOAT_BUFFER(addr,base) if(addr) {addr = (AT_FLOAT *)((int64_t)addr + (int64_t)base);}\n");
  // Add implementation specific function to the code
  addSpecificFunctions(code);

}
