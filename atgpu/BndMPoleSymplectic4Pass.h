#ifndef AT_GPU_BNDMPOLESYMPLECTIC4PASS_H
#define AT_GPU_BNDMPOLESYMPLECTIC4PASS_H
#include "StrMPoleSymplectic4Pass.h"

class BndMPoleSymplectic4Pass: public StrMPoleSymplectic4Pass {

public:
  // Construct a multipole bend pass
  explicit BndMPoleSymplectic4Pass() noexcept;
  ~BndMPoleSymplectic4Pass() noexcept override;

  // Retrieve parameters from upper layer (Python, Matlab)
  void getParameters(AbstractInterface *param, PassMethodInfo *info) override;
  void fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t *offset) override;

  // Generic code generation
  static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
  static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

  static void generateBendFringeEnter(std::string& code, PassMethodInfo *info) noexcept;
  static void generateBendFringeExit(std::string& code, PassMethodInfo *info) noexcept;

private:

  bool isBending();
  bool isDQ();

};

#endif //AT_GPU_BNDMPOLESYMPLECTIC4PASS_H
