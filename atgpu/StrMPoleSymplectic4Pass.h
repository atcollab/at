#ifndef AT_GPU_STRMPOLESYMPLECTIC4PASS_H
#define AT_GPU_STRMPOLESYMPLECTIC4PASS_H
#include "IdentityPass.h"

class StrMPoleSymplectic4Pass: public IdentityPass {

  // Construct a drift pass
  StrMPoleSymplectic4Pass() noexcept;

  // Retrieve parameters from upper layer (Python, Matlab)
  virtual void getParameters(AbstractInterface *param, PASSMETHOD_INFO *info);

  void getPolynom(AT_FLOAT *dest,AbstractInterface *param,const std::string& name,int order);

  // GPU code generation
  static void generateGPUKernel(std::string& code, PASSMETHOD_INFO *info) noexcept;
  static void generateCall(std::string& code) noexcept;

};

#endif //AT_GPU_STRMPOLESYMPLECTIC4PASS_H
