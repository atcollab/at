#ifndef AT_GPU_DRIFTPASS_H
#define AT_GPU_DRIFTPASS_H
#include "IdentityPass.h"

class DriftPass: public IdentityPass {

public:
  // Construct a drift pass
  DriftPass() noexcept;

  // Retrieve parameters from upper layer (Python, Matlab)
  void getParameters(AbstractInterface *param, PassMethodInfo *info) override;

  // Generic code generation
  static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
  static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

};

#endif //AT_GPU_DRIFTPASS_H