#ifndef AT_GPU_PASSMETHODFACTORY_H
#define AT_GPU_PASSMETHODFACTORY_H
#include "SymplecticIntegrator.h"
#include "AbstractElement.h"

class AbstractElement;
typedef AbstractElement *(*Constructor)();
typedef void (*Generator)(std::string& code,PassMethodInfo *info,SymplecticIntegrator& integrator);
typedef void (*UGenerator)(std::string& code,PassMethodInfo *info);

class PassMethodInfo {
public:
  PassMethodInfo();
  PassMethodInfo(const std::string& name,Constructor c,Generator g,UGenerator ug);
  std::string name;       // Pass method name
  Constructor create;     // Pass method constructor
  Generator generate;     // Pass method code generator
  UGenerator ugenerate;   // Pass method util functions
  bool used;              // Is pass method used ?
  bool doR1;              // Pass method use R1
  bool doR2;              // Pass method use R2
  bool doT1;              // Pass method use T1
  bool doT2;              // Pass method use T2
  bool doEAperture;       // Pass method use elliptical aperture check
  bool doRAperture;       // Pass method use rectangular aperture check
  bool doQuadEnter;       // Pass method use Quad fringe at entrance
  bool doQuadExit;        // Pass method use Quad fringe at exit
};

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

  // Generate utils code
  void generateUtilsFunctions(std::string& code);

  // Polynomial evalulation loop
  static std::string polyLoop;

private:

  SymplecticIntegrator& integrator;
  PassMethodInfo passMethodInfos[NB_PASSMETHOD_TYPE]; // Flags for pass method code generation
  std::string callCode;                               // Switch/case code for pass methods

};


#endif //AT_GPU_PASSMETHODFACTORY_H
