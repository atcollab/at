#ifndef AT_GPU_BNDMPOLESYMPLECTIC4PASS_H
#define AT_GPU_BNDMPOLESYMPLECTIC4PASS_H
#include "StrMPoleSymplectic4Pass.h"

class BndMPoleSymplectic4Pass: public StrMPoleSymplectic4Pass {

public:
  // Construct a BendMPole pass
  BndMPoleSymplectic4Pass(SymplecticIntegrator& integrator) noexcept;
  ~BndMPoleSymplectic4Pass() noexcept;

  // Retrieve parameters from upper layer (Python, Matlab)
  virtual void getParameters(AbstractInterface *param, PASSMETHOD_INFO *info);

  // GPU code generation
  static void generateGPUKernel(std::string& code, PASSMETHOD_INFO *info,SymplecticIntegrator& integrator) noexcept;
  static void generateCall(std::string& code) noexcept;

private:

  static void generateBendFringeEnter(std::string& code, PASSMETHOD_INFO *info) noexcept;
  static void generateBendFringeExit(std::string& code, PASSMETHOD_INFO *info) noexcept;

};


#endif //AT_GPU_BNDMPOLESYMPLECTIC4PASS_H
