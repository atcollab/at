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

// Abstract class to handle GPU access
class AbstractGPU {

public:

  // Return list of GPU device present on the system
  virtual std::vector<GPU_INFO> getDeviceList() = 0;

  // Return device function qualifier
  virtual void getDeviceFunctionQualifier(std::string& ftype) = 0;

  // Copy from host to device
  virtual void hostToDevice(void *dest,void *src,size_t size) = 0;

  // Allocate device memory
  virtual void allocDevice(void **dest,size_t size) = 0;

  // Compile and run the kernel
  virtual void run(std::string& code) = 0;

  // Add math function to the code
  void addUtilsFunctions(std::string& code);

  // Return handle to singleton class
  static AbstractGPU *getInstance();

  // Output code with line number
  static void outputCode(std::string& code);

private:
  static void split(std::vector<std::string> &tokens, const std::string &text, char sep);

  static AbstractGPU *gpuHandler;

};


#endif //AT_GPU_GPUWRAPPER_H
