#ifndef AT_GPU_PASSMETHODFACTORY_H
#define AT_GPU_PASSMETHODFACTORY_H
#include "AbstractElement.h"

// Singleton class to handle pass method code generation
class PassMethodFactory {

public:

  // Reset the factory
  void reset();

  // Create an element with the specified pass method
  AbstractElement *createElement(std::string& passMethod);

  // Generate pass method code
  void generatePassMethods(std::string& code);

  // Generate pass method call code
  void generatePassMethodsCalls(std::string& code);

  // Return instance of the PassMethodFactory
  static PassMethodFactory *getInstance();

private:
  PASSMETHOD_INFO passMethodInfos[NB_PASSMETHOD_TYPE];   // Flags for pass method code generation
  std::string callCode;                              // Switch/case code for pass methods
  static PassMethodFactory *handler;

};


#endif //AT_GPU_PASSMETHODFACTORY_H
