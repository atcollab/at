#ifndef AT_GPU_STRMPOLESYMPLECTIC4PASS_H
#define AT_GPU_STRMPOLESYMPLECTIC4PASS_H
#include "IdentityPass.h"
#include "SymplecticIntegrator.h"

class StrMPoleSymplectic4Pass: public IdentityPass {

public:
  // Construct a MPole pass
  explicit StrMPoleSymplectic4Pass(SymplecticIntegrator& integrator) noexcept;
  ~StrMPoleSymplectic4Pass() noexcept override;

  // Retrieve parameters from upper layer (Python, Matlab)
  void getParameters(AbstractInterface *param, PASSMETHOD_INFO *info) override;
  uint64_t getMemorySize() override;
  void fillGPUMemory(void *elemMem,void *privateMem,void *gpuMem) override;

  // GPU code generation
  static void generateGPUKernel(std::string& code, PASSMETHOD_INFO *info,SymplecticIntegrator& integrator) noexcept;
  static void generateCall(std::string& code) noexcept;
  static void generateQuadFringeEnter(std::string& code, PASSMETHOD_INFO *info) noexcept;
  static void generateQuadFringeExit(std::string& code, PASSMETHOD_INFO *info) noexcept;

protected:

  SymplecticIntegrator& integrator;
  AT_FLOAT *PolynomA;  // PolynomA
  AT_FLOAT *PolynomB;  // PolynomB
  AT_FLOAT *KickAngle; // KickAngle

private:

  bool isDrift();
  bool isQuadrupole();
  bool isSextupole();
  bool isOctupole();

};

#endif //AT_GPU_STRMPOLESYMPLECTIC4PASS_H
