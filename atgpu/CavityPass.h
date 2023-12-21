#ifndef AT_GPU_CAVITYPASS_H
#define AT_GPU_CAVITYPASS_H
#include "IdentityPass.h"

class CavityPass: public IdentityPass {

public:
  // Construct a Cavity pass
  explicit CavityPass() noexcept;
  ~CavityPass() noexcept override;

  // Retrieve parameters from upper layer (Python, Matlab)
  void getParameters(AbstractInterface *param, PassMethodInfo *info) override;

  // Generic code generation
  static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
  static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

protected:

  AT_FLOAT Frequency;

};
#endif //AT_GPU_CAVITYPASS_H
