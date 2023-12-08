#ifndef AT_GPU_ABSTRACTELEMENT_H
#define AT_GPU_ABSTRACTELEMENT_H
#include <string>
#define READER(x) x
#include "Element.h"
#include "AbstractInterface.h"

typedef struct {
  bool used;          // Is pass method used ?
  bool doR1;          // Pass method use R1 flags
  bool doR2;          // Pass method use R2 flags
  bool doT1;          // Pass method use T1 flags
  bool doT2;          // Pass method use T2 flags
  bool doEAperture;   // Pass method use Elliptical Aperture
  bool doRAperture;   // Pass method use Rectangular Aperture
} PASSMETHOD_INFO;

class AbstractInterface;

// Lattice element abstract class
class AbstractElement {

public:

  // Retrieve parameters from upper layer (Python, Matlab)
  virtual void getParameters(AbstractInterface *param,PASSMETHOD_INFO *info) = 0;

};


#endif //AT_GPU_ABSTRACTELEMENT_H
