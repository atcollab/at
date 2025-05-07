#ifndef AT_GPU_STRMPOLESYMPLECTIC4PASS_H
#define AT_GPU_STRMPOLESYMPLECTIC4PASS_H
#include "IdentityPass.h"
#include "SymplecticIntegrator.h"

class StrMPoleSymplectic4Pass: public IdentityPass {

public:
    // Construct a straight multipole pass
  explicit StrMPoleSymplectic4Pass() noexcept;
  ~StrMPoleSymplectic4Pass() noexcept override;

  // Retrieve parameters from upper layer (Python, Matlab)
  void getParameters(AbstractInterface *param, PassMethodInfo *info) override;
  uint64_t getMemorySize() override;
  void fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t *offset) override;

  // Generic code generation
  static void generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept;
  static void generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept;

  static void generateQuadFringeEnter(std::string& code, PassMethodInfo *info) noexcept;
  static void generateQuadFringeExit(std::string& code, PassMethodInfo *info) noexcept;

protected:

  AT_FLOAT *PolynomA;  // PolynomA
  AT_FLOAT *PolynomB;  // PolynomB
  AT_FLOAT *KickAngle; // KickAngle

private:

  bool isQuadrupole();
  bool isSextupole();
  bool isOctupole();

};

#endif //AT_GPU_STRMPOLESYMPLECTIC4PASS_H
