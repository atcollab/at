#include "ThinMPolePass.h"
#include "AbstractGPU.h"
#include "PassMethodFactory.h"

using namespace std;

ThinMPolePass::ThinMPolePass() noexcept : StrMPoleSymplectic4Pass() {
  BendingAngle = nullptr;
}

ThinMPolePass::~ThinMPolePass() noexcept {
}

void ThinMPolePass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  StrMPoleSymplectic4Pass::getParameters(param,info);

  elemData.Type = THINMPOLEPASS;
  param->getOptional1DArray(&BendingAngle,"BendingAngle",&BendingAngleSize);

}

void ThinMPolePass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateEnter(code,info);
  generateApertures(code,info);

  code.append("  switch(elem->SubType) {\n");
  code.append("  case 1:\n");
  code.append("    quadthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,(AT_FLOAT)1);\n");
  code.append("    break;\n");
  code.append("  case 2:\n");
  code.append("    sextuthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,(AT_FLOAT)1);\n");
  code.append("    break;\n");
  code.append("  case 3:\n");
  code.append("    octuthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,(AT_FLOAT)1);\n");
  code.append("    break;\n");
  code.append("  default:\n");
  code.append("    strthinkick(r6,elem->mpole.PolynomA,elem->mpole.PolynomB,(AT_FLOAT)1,elem->mpole.MaxOrder);\n");
  code.append("  }\n");

  code.append("  r6[1] += elem->mpole.thin.bax*r6[4];\n");
  code.append("  r6[3] -= elem->mpole.thin.bay*r6[4];\n");
  // Path lengthening
  code.append("  r6[5] -= elem->mpole.thin.bax*r6[0]-elem->mpole.thin.bay*r6[2];\n");

  generateApertures(code,info);
  generateExit(code,info);

}

void ThinMPolePass::fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t *offset) {

  // Store an offset from the beginning of gpuRing memory in ELEMENT
  // for mapping buffers in GPU memory address space (see Lattice::mapBuffers)

  StrMPoleSymplectic4Pass::fillGPUMemory(startAdd,elemMem,offset);

  if(BendingAngle) {
    if(BendingAngleSize>=1) elemData.mpole.thin.bax = BendingAngle[0];
    if(BendingAngleSize>=2) elemData.mpole.thin.bay = BendingAngle[1];
  }

  // Update modified field in buffer
  elemMem->mpole.thin.bax = elemData.mpole.thin.bax;
  elemMem->mpole.thin.bay = elemData.mpole.thin.bay;

}

void ThinMPolePass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

}