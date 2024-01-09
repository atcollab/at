#include "StrMPoleSymplectic4Pass.h"
#include "PassMethodFactory.h"
#include "AbstractGPU.h"
#include "PassMethodFactory.h"
#include <string.h>
#include <math.h>

using namespace std;

StrMPoleSymplectic4Pass::StrMPoleSymplectic4Pass() noexcept : IdentityPass() {
  PolynomA = nullptr;
  PolynomB = nullptr;
  KickAngle = nullptr;
}

StrMPoleSymplectic4Pass::~StrMPoleSymplectic4Pass() noexcept {
  delete[] PolynomA;
  delete[] PolynomB;
  delete[] KickAngle;
}

void StrMPoleSymplectic4Pass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  IdentityPass::getParameters(param,info);

  elemData.Type = STRMPOLESYMPLECTIC4PASS;
  elemData.Length = param->getDouble("Length");
  elemData.NumIntSteps = param->getOptionalInt("NumIntSteps",10);
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
    elemData.Type = DRIFTPASS;
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

void StrMPoleSymplectic4Pass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  code.append(
          "  AT_FLOAT p_norm = 1.0 / (1.0 + r6[4]);\n"
  );

  generateEnter(code,info);
  generateApertures(code,info);
  generateQuadFringeEnter(code,info);

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

}

void StrMPoleSymplectic4Pass::generateQuadFringeEnter(std::string& code, PassMethodInfo *info) noexcept {

  if(info->doQuadEnter)
    code.append(
          "  if (elem->FringeQuadEntrance) {\n"
          "    quad_fringe(r6,elem->PolynomB[1],1.0,p_norm);\n"
          "  }\n"
  );

}

void StrMPoleSymplectic4Pass::generateQuadFringeExit(std::string& code, PassMethodInfo *info) noexcept {

  if(info->doQuadExit)
    code.append(
          "  if(elem->FringeQuadExit) {\n"
          "    quad_fringe(r6,elem->PolynomB[1],-1.0,p_norm);\n"
          "  }\n"
  );

}

void StrMPoleSymplectic4Pass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  getGPUFunctionQualifier(ftype);

  // Quad
  code.append(
          ftype +
          "void quadthinkick(AT_FLOAT* r,AT_FLOAT A0,AT_FLOAT B0,AT_FLOAT K,AT_FLOAT L) {\n"
          "  r[1] -= L * (K * r[0] + B0);\n"
          "  r[3] += L * (K * r[2] + A0);\n"
          "}\n"
  );

  // Sextu (no SQ component)
  code.append(
          ftype +
          "void sextuthinkick(AT_FLOAT* r,AT_FLOAT A0,AT_FLOAT B0,AT_FLOAT K,AT_FLOAT L) {\n"
          "  r[1] -= L * (K * (r[0]*r[0]-r[2]*r[2]) + B0);\n"
          "  r[3] += L * (K * (2.0* r[0] * r[2]) + A0);\n"
          "}\n"
  );

  // Octu
  code.append(
          ftype +
          "void octuthinkick(AT_FLOAT* r,AT_FLOAT A0,AT_FLOAT B0,AT_FLOAT K,AT_FLOAT L) {\n"
          "  AT_FLOAT x2 = r[0]*r[0];\n"
          "  AT_FLOAT y2 = r[2]*r[2];\n"
          "  r[1] -= L * ((K * r[0] * (x2 - 3.0*y2)) + B0);\n"
          "  r[3] += L * ((K * r[2] * (3.0*x2 - y2)) + A0);\n"
          "}\n"
  );


  // Generic kick in straight element
  code.append(
          ftype +
          "void strthinkick(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order) {\n"
          + PassMethodFactory::polyLoop +
          "  r[1] -= L * ReSum;\n"
          "  r[3] += L * ImSum;\n"
          "}\n"
  );

  //Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam Dynamics..." by E. Forest
  code.append(
          ftype +
          "void quad_fringe(AT_FLOAT* r, AT_FLOAT b2, AT_FLOAT sign, AT_FLOAT p_norm) {\n"
          "  AT_FLOAT u = p_norm * b2 / 12.0;\n"
          "  AT_FLOAT x2 = r[0] * r[0];\n"
          "  AT_FLOAT z2 = r[2] * r[2];\n"
          "  AT_FLOAT xz = r[0] * r[2];\n"
          "  AT_FLOAT gx = u * (x2 + 3 * z2) * r[0];\n"
          "  AT_FLOAT gz = u * (z2 + 3 * x2) * r[2];\n"
          "  AT_FLOAT r1tmp = 0;\n"
          "  AT_FLOAT r3tmp = 0;\n"
          "  r[0] += sign*gx;\n"
          "  r1tmp = 3 * u * (2 * xz * r[3] - (x2 + z2) * r[1]);\n"
          "  r[2] -= sign*gz;\n"
          "  r3tmp = 3 * u * (2 * xz * r[1] - (x2 + z2) * r[3]);\n"
          "  r[5] -= sign * (gz * r[3] - gx * r[1]) * p_norm;\n"
          "  r[1] += sign*r1tmp;\n"
          "  r[3] -= sign*r3tmp;\n"
          "}\n"
  );

}
