#ifndef AT_GPU_APERTUREPASS_H
#define AT_GPU_APERTUREPASS_H
#include "AbstractElement.h"
#include "SymplecticIntegrator.h"

class AperturePass: public AbstractElement {

public:

  // Construct an Aperture pass
  AperturePass() noexcept;
  ~AperturePass() noexcept override;

  // Retrieve parameters from upper layer (Python, Matlab)
  void getParameters(AbstractInterface *param, PassMethodInfo *info) override;
  uint64_t getMemorySize() override;
  void fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t* offset) override;
  AT_FLOAT getLength() override;
  uint32_t getType() override;

  // Generic code generation
  static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
  static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

  // Code generation
  static void generateApertures(std::string& code, PassMethodInfo *info) noexcept;
  static void generateRAperture(std::string& code) noexcept;
  
protected:

  ELEMENT elemData;

private:

  AT_FLOAT *RApertures; // Rectangular transverse aperture check (xmin,xmax,ymin,ymax)

};


#endif //AT_GPU_APERTUREPASS_H
