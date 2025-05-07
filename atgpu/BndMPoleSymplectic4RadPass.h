#ifndef AT_GPU_BNDMPOLESYMPLECTIC4RADPASS_H
#define AT_GPU_BNDMPOLESYMPLECTIC4RADPASS_H

#include "BndMPoleSymplectic4Pass.h"

class BndMPoleSymplectic4RadPass: public BndMPoleSymplectic4Pass {

public:
  // Construct a multipole bend pass with radiation
  explicit BndMPoleSymplectic4RadPass() noexcept;
  ~BndMPoleSymplectic4RadPass() noexcept override;

  // Retrieve parameters from upper layer (Python, Matlab)
  void getParameters(AbstractInterface *param, PassMethodInfo *info) override;

  // Generic code generation
  static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
  static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

};

#endif //AT_GPU_BNDMPOLESYMPLECTIC4RADPASS_H
