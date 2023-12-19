#ifndef AT_GPU_ABSTRACTELEMENT_H
#define AT_GPU_ABSTRACTELEMENT_H
#include <string>
#include "Element.h"
#include "AbstractInterface.h"
#include "AbstractGPU.h"

typedef struct {
  bool used;          // Is pass method used ?
  bool doR1;          // Pass method use R1
  bool doR2;          // Pass method use R2
  bool doT1;          // Pass method use T1
  bool doT2;          // Pass method use T2
  bool doEAperture;   // Pass method use elliptical aperture check
  bool doRAperture;   // Pass method use rectangular aperture check
  bool doQuadEnter;   // Pass method use Quad fringe at entrance
  bool doQuadExit;    // Pass method use Quad fringe at exit
} PASSMETHOD_INFO;

class AbstractInterface;

// Lattice element abstract class
class AbstractElement {

public:

  virtual ~AbstractElement() noexcept = default;

  // Retrieve parameters from upper layer (Python, Matlab)
  virtual void getParameters(AbstractInterface *param,PASSMETHOD_INFO *info) = 0;

  // Get needed memory (ELEMENT not included)
  virtual uint64_t getMemorySize() = 0;

  // Fill device memory
  virtual void fillGPUMemory(void *elemMem,void *privateMem,void *gpuMem) = 0;

};


#endif //AT_GPU_ABSTRACTELEMENT_H
