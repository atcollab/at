#ifndef AT_GPU_DRIFTPASS_H
#define AT_GPU_DRIFTPASS_H
#include "IdentityPass.h"

class DriftPass: public IdentityPass {

public:
  // Construct a drift pass
  DriftPass() noexcept;

  // Retrieve parameters from upper layer (Python, Matlab)
  virtual void getParameters(AbstractInterface *param, PASSMETHOD_INFO *info);

  // GPU code generation
  static void generateGPUKernel(std::string& code, PASSMETHOD_INFO *info) noexcept;
  static void generateCall(std::string& code) noexcept;

};

#endif //AT_GPU_DRIFTPASS_H