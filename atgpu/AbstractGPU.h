#ifndef AT_GPU_GPUWRAPPER_H
#define AT_GPU_GPUWRAPPER_H
#include <cstdint>
#include <string>
#include <vector>

// Abstract class for GPU context
// Abstract OpenCL or CUDA API

// GPU device information

typedef struct {
  std::string name;        // GPU Name
  std::string version;     // Compute capabilities (CUDA)@guney1911
  uint32_t    mpNumber;    // Multi processor number
  std::string platform;    // Platform name
} GPU_INFO;

// MACRO to define GPU Context interface

#define GPU_CONTEXT_DECL(_funcQualifier,funcInit)                                  \
/* Empty kernel parameter */                                                       \
_funcQualifier void resetArg()funcInit;                                            \
/* Set kernel parameter */                                                         \
_funcQualifier void addArg(size_t argSize,void *value)funcInit;                    \
/* Run the kernel */                                                               \
_funcQualifier void run(uint64_t nbThread)funcInit;                                \
/* Compile and load the kernel */                                                  \
_funcQualifier void compile(std::string& code)funcInit;                            \
/* Map memory buffer */                                                            \
_funcQualifier void mapBuffer(void **ring,uint32_t nbElement)funcInit;             \
/* Copy from host to device */                                                     \
_funcQualifier void hostToDevice(void *dest,void *src,size_t size)funcInit;        \
/* Copy from device to host */                                                     \
_funcQualifier void deviceToHost(void *dest,void *src,size_t size)funcInit;        \
/* Allocate device memory */                                                       \
_funcQualifier void allocDevice(void **dest,size_t size)funcInit;                  \
/* Free device memory */                                                           \
_funcQualifier void freeDevice(void *dest)funcInit;                                \
/* Return device name */                                                           \
_funcQualifier std::string& name()funcInit;                                        \
/* Return number fo core */                                                        \
_funcQualifier int coreNumber()funcInit;

// MACRO to define GPU Context implementation

#define GPU_CONTEXT_IMPL() GPU_CONTEXT_DECL(,override)

// MACRO to define GPU singleton class interface (for context creation)

#define GPU_DECL(_funcQualifier,_funcInit)                                         \
/* Return list of GPU device present on the system */                              \
_funcQualifier std::vector<GPU_INFO> getDeviceList()_funcInit;                     \
/* Create a context for the given GPU */                                           \
_funcQualifier GPUContext *createContext(int devId)_funcInit;                      \
/* Return device function qualifier */                                             \
_funcQualifier void getDeviceFunctionQualifier(std::string& ftype)_funcInit;       \
/* Return kernel function qualifier */                                             \
_funcQualifier void getKernelFunctionQualifier(std::string& ftype)_funcInit;       \
/* Return kernel function qualifier */                                             \
_funcQualifier void getGlobalQualifier(std::string& ftype)_funcInit;               \
/* Return command to compute the thread id */                                      \
_funcQualifier void getThreadId(std::string& command)_funcInit;                    \
/* Format a double */                                                              \
_funcQualifier std::string formatFloat(double *f)_funcInit;                        \
/* Add implementation specific function to the code */                             \
_funcQualifier void addSpecificFunctions(std::string& code)_funcInit;

// MACRO to define GPU singleton class implementation

#define GPU_IMPL() GPU_DECL(,override)


class GPUContext {
public:
  virtual ~GPUContext() = default;
  // Declare GPU context interface as pure virtual function
  GPU_CONTEXT_DECL(virtual,=0)
};

// Abstract class to handle GPU access
class AbstractGPU {

public:

  GPU_DECL(virtual,=0)

  // Add util functions to the code
  void addUtilsFunctions(std::string &code);

  // Return implementation name
  std::string& implName();

  // Return handle to singleton class
  static AbstractGPU *getInstance();

  // Output code with line number
  static void outputCode(std::string& code);

  // Get number of sec since instantiation of this singleton class
  static double get_ticks();

protected:

  // Initialisation error
  std::string initErrorStr;
  std::string implementationStr;

private:

  static void initTimer();
  static AbstractGPU *gpuHandler;

};


#endif //AT_GPU_GPUWRAPPER_H
