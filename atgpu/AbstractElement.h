#ifndef AT_GPU_ABSTRACTELEMENT_H
#define AT_GPU_ABSTRACTELEMENT_H
#include <string>
#include <cstdint>
#include "Element.h"
#include "AbstractInterface.h"
#include "AbstractGPU.h"

class PassMethodInfo;
class AbstractInterface;

// Lattice element abstract class
class AbstractElement {

public:

  virtual ~AbstractElement() noexcept = default;

  // Retrieve parameters from upper layer (Python, Matlab)
  virtual void getParameters(AbstractInterface *param,PassMethodInfo *info) = 0;

  // Get needed memory (ELEMENT not included)
  virtual uint64_t getMemorySize() = 0;

  // Fill device memory
  virtual void fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t* offset) = 0;

  // Get element length
  virtual AT_FLOAT getLength() = 0;

  // Get element type
  virtual uint32_t getType() = 0;

  // Post initialisation
  virtual void postInit(RING_PARAM *param) {};

};


#endif //AT_GPU_ABSTRACTELEMENT_H
