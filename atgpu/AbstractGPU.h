#ifndef AT_GPU_GPUWRAPPER_H
#define AT_GPU_GPUWRAPPER_H

#include <string>
#include <vector>

typedef struct {
  std::string name;        // GPU Name
  std::string version;     // Compute capabilities
  uint32_t    smNumber;    // Stream processor number
  uint32_t    mpNumber;    // Multi processor number
} GPU_INFO;

// Abstract class for GPU context

class GPUContext {

public:

  virtual ~GPUContext() = default;

  // Empty kernel parameter
  virtual void resetArg() = 0;

  // Set kernel parameter
  virtual void addArg(size_t argSize,void *value) = 0;

  // Run the kernel
  virtual void run(uint32_t gridSize,uint64_t nbThread) = 0;

  // Compile and load the kernel
  virtual void compile(std::string& code) = 0;

  // Copy from host to device
  virtual void hostToDevice(void *dest,void *src,size_t size) = 0;

  // Copy from device to host
  virtual void deviceToHost(void *dest,void *src,size_t size) = 0;

  // Allocate device memory
  virtual void allocDevice(void **dest,size_t size) = 0;

  // Free device memory
  virtual void freeDevice(void *dest) = 0;

  // Return device name
  virtual std::string& name() = 0;

};

// Abstract class to handle GPU access
class AbstractGPU {

public:

  // Return list of GPU device present on the system
  virtual std::vector<GPU_INFO> getDeviceList() = 0;

  // Create a context for the given GPU
  virtual GPUContext *createContext(int devId) = 0;

  // Return device function qualifier
  virtual void getDeviceFunctionQualifier(std::string& ftype) = 0;

  // Add utils function to the code
  void addUtilsFunctions(std::string& code);

  // Return handle to singleton class
  static AbstractGPU *getInstance();

  // Output code with line number
  static void outputCode(std::string& code);

  // Get number of sec since instantiation of this singleton class
  static double get_ticks();

private:

  static void split(std::vector<std::string> &tokens, const std::string &text, char sep);
  static void initTimer();

  static AbstractGPU *gpuHandler;

};


#endif //AT_GPU_GPUWRAPPER_H
