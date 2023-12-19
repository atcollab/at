#include "StrMPoleSymplectic4Pass.h"
#include "AbstractGPU.h"
#include <iostream>
#include <string.h>
#include <math.h>

using namespace std;

StrMPoleSymplectic4Pass::StrMPoleSymplectic4Pass(SymplecticIntegrator& integrator) noexcept : IdentityPass(),
integrator(integrator) {
  PolynomA = nullptr;
  PolynomB = nullptr;
  KickAngle = nullptr;
}

StrMPoleSymplectic4Pass::~StrMPoleSymplectic4Pass() noexcept {
  delete[] PolynomA;
  delete[] PolynomB;
  delete[] KickAngle;
}

void StrMPoleSymplectic4Pass::getParameters(AbstractInterface *param, PASSMETHOD_INFO *info) {

  // Retrieve param from super class
  IdentityPass::getParameters(param,info);

  elemData.Type = MPOLE;
  elemData.Length = param->getDouble("Length");
  elemData.NumIntSteps = param->getInt("NumIntSteps");
  elemData.SL = elemData.Length / (AT_FLOAT)elemData.NumIntSteps;
  elemData.MaxOrder = param->getInt("MaxOrder");

  param->get1DArray(&PolynomA,"PolynomA",elemData.MaxOrder+1);
  param->get1DArray(&PolynomB,"PolynomB",elemData.MaxOrder+1);
  param->getOptional1DArray(&KickAngle,"KickAngle",2);
  if( KickAngle ) {
    PolynomB[0] -= sin(KickAngle[0]) / elemData.Length;
    PolynomA[0] += sin(KickAngle[1]) / elemData.Length;
  }

  elemData.FringeQuadEntrance = param->getOptionalInt("FringeQuadEntrance", 0);
  elemData.FringeQuadExit = param->getOptionalInt("FringeQuadExit", 0);
  if( elemData.MaxOrder>=1 && PolynomB[1]==0.0 ) {
    // No quad strength
    elemData.FringeQuadEntrance = false;
    elemData.FringeQuadExit = false;
  }

  if ( isDrift() ) {
    // All polynom coefficients are null
    elemData.Type = DRIFT;
  } else if( isQuadrupole() ) {
    elemData.SubType = 1;
    elemData.K = PolynomB[1];
  } else if ( isSextupole() ) {
    elemData.SubType = 2;
    elemData.K = PolynomB[2];
  } else if ( isOctupole() ) {
    elemData.SubType = 3;
    elemData.K = PolynomB[3];
  }

  info->doQuadEnter |= (elemData.FringeQuadEntrance!=0);
  info->doQuadExit |= (elemData.FringeQuadExit!=0);

}

uint64_t StrMPoleSymplectic4Pass::getMemorySize() {

  uint64_t sum = IdentityPass::getMemorySize();
  if(PolynomA) sum += (elemData.MaxOrder + 1) * sizeof(AT_FLOAT);
  if(PolynomB) sum += (elemData.MaxOrder + 1) * sizeof(AT_FLOAT);
  return sum;

}

void StrMPoleSymplectic4Pass::fillGPUMemory(void *elemMem,void *privateMem,void *gpuMem) {

  uint64_t privSize = IdentityPass::getMemorySize();
  unsigned char *privPtr = (unsigned char *)privateMem;
  AT_FLOAT *dest = (AT_FLOAT *)(privPtr + privSize);
  AT_FLOAT *destGPU = (AT_FLOAT *)((unsigned char *)gpuMem + privSize);

  if(PolynomA) {
    elemData.PolynomA = destGPU;
    memcpy(dest,PolynomA,(elemData.MaxOrder + 1)*sizeof(AT_FLOAT));
    dest += (elemData.MaxOrder + 1);
    destGPU += (elemData.MaxOrder + 1);
  }
  if(PolynomB) {
    elemData.PolynomB = destGPU;
    memcpy(dest,PolynomB,(elemData.MaxOrder + 1)*sizeof(AT_FLOAT));
    dest += (elemData.MaxOrder + 1);
    destGPU += (elemData.MaxOrder + 1);
  }

  IdentityPass::fillGPUMemory(elemMem,privateMem,gpuMem);

}

bool StrMPoleSymplectic4Pass::isDrift() {
  bool isDr = true;
  for(int i=0;i<=elemData.MaxOrder;i++)
    isDr &= PolynomA[i]==0.0 && PolynomB[i]==0.0;
  return isDr;
}

bool StrMPoleSymplectic4Pass::isQuadrupole() {
  return elemData.MaxOrder==1 && PolynomA[1]==0.0;
}

bool StrMPoleSymplectic4Pass::isSextupole() {
  return elemData.MaxOrder==2 && PolynomA[2]==0.0 &&
         PolynomA[1]==0.0 && PolynomB[1]==0.0;
}

bool StrMPoleSymplectic4Pass::isOctupole() {
  return elemData.MaxOrder==3 && PolynomA[3]==0.0 &&
         PolynomA[2]==0.0 && PolynomB[2]==0.0 &&
         PolynomA[1]==0.0 && PolynomB[1]==0.0;
}

void StrMPoleSymplectic4Pass::generateGPUKernel(std::string& code, PASSMETHOD_INFO *info,SymplecticIntegrator& integrator) noexcept {

  AbstractGPU *gpu = AbstractGPU::getInstance();
  string ftype;
  gpu->getDeviceFunctionQualifier(ftype);
  if(!ftype.empty()) ftype.append(" ");

  code.append( ftype + "void StrMPoleSymplectic4Pass(AT_FLOAT* r6,ELEMENT* elem) {\n");
  code.append(
          "  AT_FLOAT p_norm = 1.0 / (1.0 + r6[4]);\n"
  );

  generateEnter(code,info);
  generateApertures(code,info);
  generateQuadFringeEnter(code,info);

  // Kick/Drift methods are defined in PassMethodFactory
  integrator.resetMethods();
  // Default straight magnet
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("strthinkick(r6,elem->PolynomA,elem->PolynomB,%STEP%,elem->MaxOrder)");
  // Pure quad
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("quadthinkick(r6,elem->PolynomA[0],elem->PolynomB[0],elem->K,%STEP%)");
  // Pure sextu
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("sextuthinkick(r6,elem->PolynomA[0],elem->PolynomB[0],elem->K,%STEP%)");
  // Pure octu
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("octuthinkick(r6,elem->PolynomA[0],elem->PolynomB[0],elem->K,%STEP%)");

  integrator.generateCode(code);

  generateQuadFringeExit(code,info);
  generateApertures(code,info);
  generateExit(code,info);
  code.append("}\n");

}


void StrMPoleSymplectic4Pass::generateCall(std::string& code) noexcept {

  code.append(
          "      case MPOLE:\n"
          "        StrMPoleSymplectic4Pass(r6,elemPtr);\n"
          "        break;\n"
  );

}

void StrMPoleSymplectic4Pass::generateQuadFringeEnter(std::string& code, PASSMETHOD_INFO *info) noexcept {

  if(info->doQuadEnter)
    code.append(
          "  if (elem->FringeQuadEntrance) {\n"
          "    quad_fringe(r6,elem->PolynomB[1],1.0,p_norm);\n"
          "  }\n"
  );

}

void StrMPoleSymplectic4Pass::generateQuadFringeExit(std::string& code, PASSMETHOD_INFO *info) noexcept {

  if(info->doQuadExit)
    code.append(
          "  if(elem->FringeQuadExit) {\n"
          "    quad_fringe(r6,elem->PolynomB[1],-1.0,p_norm);\n"
          "  }\n"
  );

}
