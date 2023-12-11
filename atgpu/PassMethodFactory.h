#ifndef AT_GPU_PASSMETHODFACTORY_H
#define AT_GPU_PASSMETHODFACTORY_H
#include "SymplecticIntegrator.h"
#include "AbstractElement.h"

class AbstractElement;

// Class to handle pass method code generation
class PassMethodFactory {

public:

  // Construct a factory
  explicit PassMethodFactory(SymplecticIntegrator& integrator) noexcept;

  // Create an element with the specified pass method
  AbstractElement *createElement(std::string& passMethod);

  // Generate pass method code
  void generatePassMethods(std::string& code);

  // Generate pass method call code
  void generatePassMethodsCalls(std::string& code);

private:

  SymplecticIntegrator& integrator;
  PASSMETHOD_INFO passMethodInfos[NB_PASSMETHOD_TYPE];   // Flags for pass method code generation
  std::string callCode;                                  // Switch/case code for pass methods

};


#endif //AT_GPU_PASSMETHODFACTORY_H
