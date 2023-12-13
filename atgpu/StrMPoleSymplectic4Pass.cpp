#include "StrMPoleSymplectic4Pass.h"
#include "AbstractGPU.h"
#include <iostream>

using namespace std;

StrMPoleSymplectic4Pass::StrMPoleSymplectic4Pass(SymplecticIntegrator& integrator) noexcept : IdentityPass(),
integrator(integrator) {
  NormD = new AT_FLOAT[integrator.nbCoefficients];
  NormK = new AT_FLOAT[integrator.nbCoefficients];
  PolynomA = nullptr;
  PolynomB = nullptr;
  KickAngle = nullptr;
}

StrMPoleSymplectic4Pass::~StrMPoleSymplectic4Pass() noexcept {
  delete[] NormD;
  delete[] NormK;
  delete[] PolynomA;
  delete[] PolynomB;
  delete[] KickAngle;
}

void StrMPoleSymplectic4Pass::getParameters(AbstractInterface *param, PASSMETHOD_INFO *info) {

  // Retrieve param from super class
  IdentityPass::getParameters(param,info);

  elemData.Type = MPOLE;
  elemData.Length = param->getDouble("Length");
  elemData.SL = elemData.Length / (AT_FLOAT)elemData.NumIntSteps;
  elemData.MaxOrder = param->getInt("MaxOrder");

  int nbCoef = integrator.nbCoefficients;
  for(int i=0;i<nbCoef;i++) {
    NormD[i] = elemData.SL * integrator.c[i];
    NormK[i] = elemData.SL * integrator.d[i];
  }

  param->get1DArray(&PolynomA,"PolynomA",elemData.MaxOrder+1);
  param->get1DArray(&PolynomB,"PolynomB",elemData.MaxOrder+1);
  param->getOptional1DArray(&KickAngle,"KickAngle",2);

  elemData.FringeQuadEntrance = param->getOptionalInt("FringeQuadEntrance", 0);
  elemData.FringeQuadExit = param->getOptionalInt("FringeQuadExit", 0);
  if( elemData.MaxOrder>=1 && PolynomB[1]==0.0 ) {
    // No quad strength
    elemData.FringeQuadEntrance = false;
    elemData.FringeQuadExit = false;
  }

  if( isQuadrupole() )
    elemData.SubType = 1;
  else if ( isSextupole() )
    elemData.SubType = 2;
  else if ( isOctupole() )
    elemData.SubType = 3;

  info->doQuadEnter |= (elemData.FringeQuadEntrance!=0);
  info->doQuadExit |= (elemData.FringeQuadExit!=0);
  info->doKickAngle |= (elemData.KickAngle != nullptr);

}

uint64_t StrMPoleSymplectic4Pass::getMemorySize() {

  uint64_t sum = IdentityPass::getMemorySize();
  if(PolynomA) sum += (elemData.MaxOrder + 1) * sizeof(AT_FLOAT);
  if(PolynomB) sum += (elemData.MaxOrder + 1) * sizeof(AT_FLOAT);
  if(KickAngle) sum += 2 * sizeof(AT_FLOAT);
  sum += integrator.nbCoefficients * 2 * sizeof(AT_FLOAT);
  return sum;

}

void StrMPoleSymplectic4Pass::fillGPUMemory(GPUContext *gpu,void *elemMem,void *privateMem) {

  uint64_t privSize = IdentityPass::getMemorySize();
  unsigned char *privPtr = (unsigned char *)privateMem;
  AT_FLOAT *dest = (AT_FLOAT *)(privPtr + privSize);

  if(PolynomA) {
    elemData.PolynomA = dest;
    gpu->hostToDevice(dest,PolynomA,(elemData.MaxOrder + 1)*sizeof(AT_FLOAT));
    dest += (elemData.MaxOrder + 1);
  }
  if(PolynomB) {
    elemData.PolynomB = dest;
    gpu->hostToDevice(dest,PolynomB,(elemData.MaxOrder + 1)*sizeof(AT_FLOAT));
    dest += (elemData.MaxOrder + 1);
  }
  if(KickAngle) {
    elemData.KickAngle = dest;
    gpu->hostToDevice(dest,KickAngle,2*sizeof(AT_FLOAT));
    dest += 2;
  }

  elemData.NormD = dest;
  gpu->hostToDevice(dest,NormD,integrator.nbCoefficients * sizeof(AT_FLOAT));
  elemData.NormK = dest + integrator.nbCoefficients;
  gpu->hostToDevice(dest,NormK,integrator.nbCoefficients * sizeof(AT_FLOAT));

  IdentityPass::fillGPUMemory(gpu,elemMem,privateMem);

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
  generateKickAngle(code,info);
  generateQuadFringeEnter(code,info);

  // Generate switch/case for subtype (pure Quad,Sextu,Octu)
  code.append("  switch(elem->SubType) {\n");
  for(int subType=0;subType<4;subType++)
    generateIntegrator(code,subType,info,integrator);
  code.append("  }\n");

  generateQuadFringeExit(code,info);
  generateKickAngleRestore(code,info);
  generateApertures(code,info);
  generateExit(code,info);
  code.append("}\n");

}

void StrMPoleSymplectic4Pass::generateIntegrator(std::string& code, int subType, PASSMETHOD_INFO *info,SymplecticIntegrator& integrator) noexcept {

  // Integrator loop
  string sTStr = to_string(subType);
  code.append("  case " + sTStr + ":\n");
  integrator.generateCode(code,"elem->NumIntSteps","fastdrift","strthinkick" + ((subType>0)?to_string(subType):""),"");
  code.append("    break;\n");

}

void StrMPoleSymplectic4Pass::generateCall(std::string& code) noexcept {

  code.append(
          "      case MPOLE:\n"
          "        StrMPoleSymplectic4Pass(r6,elemPtr);\n"
          "        break;\n"
  );

}

void StrMPoleSymplectic4Pass::generateKickAngle(std::string& code, PASSMETHOD_INFO *info) noexcept {

  if(info->doKickAngle)
  code.append(
          "  AT_FLOAT A0 = elem->PolynomA[0];\n"
          "  AT_FLOAT B0 = elem->PolynomB[0];\n"
          "  if (elem->doKickAngle) {\n"
          "    elem->PolynomB[0] -= sin(elem->KickAngle[0])/elemData.SL;\n"
          "    elem->PolynomA[0] += sin(elem->KickAngle[1])/elemData.SL;\n"
          "  }\n"
  );

}

void StrMPoleSymplectic4Pass::generateKickAngleRestore(std::string& code, PASSMETHOD_INFO *info) noexcept {

  if(info->doKickAngle)
  code.append(
          "  if (elem->doKickAngle) {\n"
          "    elem->PolynomB[0] = B0;\n"
          "  }\n"
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
