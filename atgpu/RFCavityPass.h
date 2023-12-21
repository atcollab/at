#ifndef AT_GPU_RFCAVITYPASS_H
#define AT_GPU_RFCAVITYPASS_H
#include "CavityPass.h"

class RFCavityPass: public CavityPass {

public:
  // Construct a Cavity pass
  explicit RFCavityPass() noexcept;
  ~RFCavityPass() noexcept override;

  // Retrieve parameters from upper layer (Python, Matlab)
  void getParameters(AbstractInterface *param, PassMethodInfo *info) override;

  void postInit(RING_PARAM *param);

  // Generic code generation
  static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
  static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

};
#endif //AT_GPU_RFCAVITYPASS_H
