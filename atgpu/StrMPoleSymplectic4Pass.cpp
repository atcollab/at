#include "StrMPoleSymplectic4Pass.h"
#include "AbstractGPU.h"

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

  /*
  elemData.FringeBendEntrance = param->getOptionalInt("FringeBendEntrance",0);
  if(elemData.FringeBendEntrance) {
    elemData.EntranceAngle = param->getOptionalDouble("EntranceAngle", 0);
    elemData.tgEntranceAngle = tan(elemData.EntranceAngle);
  }
  if(elemData.FringeBendExit) {
    elemData.ExitAngle = param->getOptionalDouble("EntranceAngle", 0);
    elemData.tgEntranceAngle = tan(elemData.EntranceAngle);
  }
  */

  elemData.FringeQuadEntrance = param->getOptionalInt("FringeQuadEntrance", 0);
  elemData.FringeQuadExit = param->getOptionalInt("FringeQuadExit", 0);
  if( elemData.MaxOrder>=1 && PolynomB[1]==0.0 ) {
    // No quad strength
    elemData.FringeQuadEntrance = false;
    elemData.FringeQuadExit = false;
  }

  info->doQuadEnter |= (elemData.FringeQuadEntrance!=0);
  info->doQuadExit |= (elemData.FringeQuadExit!=0);
  info->doKickAngle |= (elemData.KickAngle != nullptr);

}

uint64_t StrMPoleSymplectic4Pass::getMemorySize() {

  uint64_t sum = 0;
  if(PolynomA) sum += (elemData.MaxOrder + 1) * sizeof(AT_FLOAT);
  if(PolynomB) sum += (elemData.MaxOrder + 1) * sizeof(AT_FLOAT);
  if(KickAngle) sum += 2 * sizeof(AT_FLOAT);
  sum += integrator.nbCoefficients * 2 * sizeof(AT_FLOAT);
  return sum;

}

void StrMPoleSymplectic4Pass::fillGPUMemory(void *deviceMem) {

  AbstractGPU *gpu = AbstractGPU::getInstance();
  AT_FLOAT *dest = (AT_FLOAT *)deviceMem;

  if(PolynomA) {
    elemData.PolynomA = dest;
    gpu->hostToDevice(dest,PolynomA,6*6*sizeof(AT_FLOAT));
    dest += (elemData.MaxOrder + 1);
  }
  if(PolynomB) {
    elemData.PolynomB = dest;
    gpu->hostToDevice(dest,PolynomB,6*6*sizeof(AT_FLOAT));
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

}

void StrMPoleSymplectic4Pass::generateGPUKernel(std::string& code, PASSMETHOD_INFO *info,SymplecticIntegrator& integrator) noexcept {

  AbstractGPU *gpu = AbstractGPU::getInstance();
  string ftype;
  gpu->getDeviceFunctionQualifier(ftype);
  if(!ftype.empty()) ftype.append(" ");

  code.append( ftype + "void StrMPoleSymplectic4Pass(AT_FLOAT* r6,ELEMENT* elem) {\n");
  generateEnter(code,info);
  generateApertures(code,info);

  code.append(
          "  AT_FLOAT p_norm = 1.0 / (1.0 + r6[4]);\n"
  );

  if(info->doKickAngle) generateKickAngle(code,info);
  if(info->doQuadEnter) generateQuadFringeEnter(code,info);


  // Integrator loop
  code.append(
          "  for(int m = 0; m < elem->NumIntSteps; m++) {\n"
  );

  for(int i=0;i<integrator.nbCoefficients;i++) {

    if( integrator.c[i]!=0.0 )
      code.append("    fastdrift(r6,elem->NormD["+to_string(i)+"] * p_norm,p_norm);\n");
    if( integrator.d[i]!=0.0 )
      code.append("    strthinkick(r6,elem->PolynomA,elem->PolynomB,elem->NormK["+to_string(i)+"],elem->MaxOrder);\n");

  }

  code.append("  }\n");

   generateQuadFringeExit(code,info);
   generateKickAngleRestore(code,info);
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
