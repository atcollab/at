#ifndef AT_GPU_STRMPOLESYMPLECTIC4RADPASS_H
#define AT_GPU_STRMPOLESYMPLECTIC4RADPASS_H

#include "StrMPoleSymplectic4Pass.h"

class StrMPoleSymplectic4RadPass: public StrMPoleSymplectic4Pass {

public:
  // Construct a straight multipole pass with radiation
  explicit StrMPoleSymplectic4RadPass() noexcept;
  ~StrMPoleSymplectic4RadPass() noexcept override;

  // Retrieve parameters from upper layer (Python, Matlab)
  void getParameters(AbstractInterface *param, PassMethodInfo *info) override;

  // Generic code generation
  static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
  static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

};

#endif //AT_GPU_STRMPOLESYMPLECTIC4RADPASS_H
