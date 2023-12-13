#ifndef AT_GPU_SYMPLECTICINTEGRATOR_H
#define AT_GPU_SYMPLECTICINTEGRATOR_H
#include <string>
#include "Element.h"

class SymplecticIntegrator {

public:

  // Construct a symplectic integrator of the given type
  SymplecticIntegrator(int type);
  ~SymplecticIntegrator();

  // Coefficients
  int nbCoefficients;
  AT_FLOAT *c; // Drift
  AT_FLOAT *d; // Kick

  void generateCode(std::string& code,const std::string& count,const std::string& driftMethod,
                    const std::string& kickMethod, const std::string& exKickParam);

private:
  void allocate(int nb);

};

#endif //AT_GPU_SYMPLECTICINTEGRATOR_H
