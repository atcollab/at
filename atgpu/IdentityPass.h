#ifndef AT_GPU_IDENTITYPASS_H
#define AT_GPU_IDENTITYPASS_H
#include "AbstractElement.h"
#include "SymplecticIntegrator.h"

class IdentityPass: public AbstractElement {

public:

  // Construct an identity pass
  IdentityPass() noexcept;
  ~IdentityPass() noexcept override;

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
  static void generateEnter(std::string& code, PassMethodInfo *info) noexcept;
  static void generateExit(std::string& code, PassMethodInfo *info) noexcept;
  static void generateApertures(std::string& code, PassMethodInfo *info) noexcept;
  static void generateEAperture(std::string& code) noexcept;
  static void generateRAperture(std::string& code) noexcept;
  static void generateR(std::string& code,const std::string& pname) noexcept;
  static void generateT(std::string& code,const std::string& pname) noexcept;
  
protected:

  ELEMENT elemData;

private:

  AT_FLOAT *R1;    // Enter 6x6 transformation matrix
  AT_FLOAT *R2;    // Exit 6x6 transformation matrix
  AT_FLOAT *T1;    // Enter 6D vector translation
  AT_FLOAT *T2;    // Exit 6D vector translation
  AT_FLOAT *EApertures; // Elliptical transverse aperture check (xradius,yradius)
  AT_FLOAT *RApertures; // Rectangular transverse aperture check (xmin,xmax,ymin,ymax)

};


#endif //AT_GPU_IDENTITYPASS_H
