#include "CorrectorPass.h"
#include <cstring>

using namespace std;

CorrectorPass::CorrectorPass() noexcept : IdentityPass() {
  KickAngle = nullptr;
}

void CorrectorPass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  IdentityPass::getParameters(param,info);

  elemData.Type = CORRECTORPASS;
  elemData.Length = param->getDouble("Length");
  param->get1DArray(&KickAngle,"KickAngle",2);

}

uint64_t CorrectorPass::getMemorySize() {

  uint64_t sum = IdentityPass::getMemorySize();
  sum += 2 * sizeof(AT_FLOAT);
  return sum;

}

void CorrectorPass::fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t *offset) {

  // Store an offset from the beginning of gpuRing memory in ELEMENT
  // for mapping buffers in GPU memory address space (see Lattice::mapBuffers)
  IdentityPass::fillGPUMemory(startAdd,elemMem,offset);

  elemData.cp.KickAngle = (AT_FLOAT *)(*offset);
  AT_FLOAT *P = (AT_FLOAT *)(startAdd+*offset);
  memcpy(P,KickAngle,2*sizeof(AT_FLOAT));
  *offset += 2*sizeof(AT_FLOAT);

  // Update modified field in buffer
  elemMem->SubType = elemData.SubType;
  elemMem->cp.KickAngle = elemData.cp.KickAngle;

}

void CorrectorPass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateEnter(code,info);
  generateApertures(code,info);
  code.append("  AT_FLOAT p_norm = PNORM(r6[4]);\n");
  code.append("  AT_FLOAT NormL = elem->Length*p_norm;\n");
  code.append("  AT_FLOAT xkick = elem->cp.KickAngle[0];\n");
  code.append("  AT_FLOAT ykick = elem->cp.KickAngle[1];\n");
  code.append("  r6[5] += NormL*p_norm*(xkick*xkick/3 + ykick*ykick/3 +\n"
              "                         r6[1]*r6[1] + r6[3]*r6[3] +\n"
              "                         r6[1]*xkick + r6[3]*ykick)/2;\n");
  code.append("  r6[0] += NormL*(r6[1]+xkick/2);\n");
  code.append("  r6[1] += xkick;\n");
  code.append("  r6[2] += NormL*(r6[3]+ykick/2);\n");
  code.append("  r6[3] += ykick;\n");
  generateApertures(code,info);
  generateExit(code,info);

}

void CorrectorPass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {
}
