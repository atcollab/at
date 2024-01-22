#include "StrMPoleSymplectic4Pass.h"
#include "PassMethodFactory.h"
#include "AbstractGPU.h"
#include <cstring>
#include <cmath>

using namespace std;

StrMPoleSymplectic4Pass::StrMPoleSymplectic4Pass() noexcept : IdentityPass() {
  PolynomA = nullptr;
  PolynomB = nullptr;
  KickAngle = nullptr;
}

StrMPoleSymplectic4Pass::~StrMPoleSymplectic4Pass() noexcept {
}

void StrMPoleSymplectic4Pass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  IdentityPass::getParameters(param,info);

  elemData.Type = STRMPOLESYMPLECTIC4PASS;
  elemData.Length = param->getDouble("Length");
  elemData.mpole.NumIntSteps = param->getOptionalInt("NumIntSteps",10);
  elemData.SL = elemData.Length / (AT_FLOAT)elemData.mpole.NumIntSteps;
  elemData.mpole.MaxOrder = param->getInt("MaxOrder");

  param->get1DArray(&PolynomA,"PolynomA",elemData.mpole.MaxOrder+1);
  param->get1DArray(&PolynomB,"PolynomB",elemData.mpole.MaxOrder+1);
  param->getOptional1DArray(&KickAngle,"KickAngle",2);

  elemData.mpole.FringeQuadEntrance = param->getOptionalInt("FringeQuadEntrance", 0);
  elemData.mpole.FringeQuadExit = param->getOptionalInt("FringeQuadExit", 0);

  info->doQuadEnter |= (elemData.mpole.FringeQuadEntrance!=0);
  info->doQuadExit |= (elemData.mpole.FringeQuadExit!=0);

}

uint64_t StrMPoleSymplectic4Pass::getMemorySize() {

  uint64_t sum = IdentityPass::getMemorySize();
  if(PolynomA) sum += (elemData.mpole.MaxOrder + 1) * sizeof(AT_FLOAT);
  if(PolynomB) sum += (elemData.mpole.MaxOrder + 1) * sizeof(AT_FLOAT);
  return sum;

}

void StrMPoleSymplectic4Pass::fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t *offset) {

  // Store an offset from the beginning of gpuRing memory in ELEMENT
  // for mapping buffers in GPU memory address space (see Lattice::mapBuffers)

  IdentityPass::fillGPUMemory(startAdd,elemMem,offset);

  if( isQuadrupole() ) {
    elemData.SubType = 1;
    elemData.mpole.K = PolynomB[1];
  } else if ( isSextupole() ) {
    elemData.SubType = 2;
    elemData.mpole.K = PolynomB[2];
  } else if ( isOctupole() ) {
    elemData.SubType = 3;
    elemData.mpole.K = PolynomB[3];
  } else {
    elemData.SubType = 0;
  }

  if(PolynomA) {
    elemData.mpole.PolynomA = (AT_FLOAT *)(*offset);
    AT_FLOAT *PA = (AT_FLOAT *)(startAdd+*offset);
    memcpy(PA,PolynomA,(elemData.mpole.MaxOrder + 1)*sizeof(AT_FLOAT));
    if( KickAngle )
      PA[0] += sin(KickAngle[1]) / elemData.Length;
    *offset += (elemData.mpole.MaxOrder + 1)*sizeof(AT_FLOAT);
  }
  if(PolynomB) {
    elemData.mpole.PolynomB = (AT_FLOAT *)(*offset);
    AT_FLOAT *PB = (AT_FLOAT *)(startAdd+*offset);
    memcpy(PB,PolynomB,(elemData.mpole.MaxOrder + 1)*sizeof(AT_FLOAT));
    if( KickAngle )
      PB[0] -= sin(KickAngle[0]) / elemData.Length;
    *offset += (elemData.mpole.MaxOrder + 1)*sizeof(AT_FLOAT);
  }

  // Update modified field in buffer
  elemMem->SubType = elemData.SubType;
  elemMem->mpole.K = elemData.mpole.K;
  elemMem->mpole.PolynomA = elemData.mpole.PolynomA;
  elemMem->mpole.PolynomB = elemData.mpole.PolynomB;

}

bool StrMPoleSymplectic4Pass::isQuadrupole() {
  return elemData.mpole.MaxOrder==1 && PolynomA[1]==0.0;
}

bool StrMPoleSymplectic4Pass::isSextupole() {
  return elemData.mpole.MaxOrder==2 && PolynomA[2]==0.0 &&
         PolynomA[1]==0.0 && PolynomB[1]==0.0;
}

bool StrMPoleSymplectic4Pass::isOctupole() {
  return elemData.mpole.MaxOrder==3 && PolynomA[3]==0.0 &&
         PolynomA[2]==0.0 && PolynomB[2]==0.0 &&
         PolynomA[1]==0.0 && PolynomB[1]==0.0;
}

void StrMPoleSymplectic4Pass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  code.append(
          "  AT_FLOAT p_norm = PNORM(r6[4]);\n"
  );

  generateEnter(code,info);
  generateApertures(code,info);
  generateQuadFringeEnter(code,info);

  integrator.resetMethods();
  // Default straight magnet
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("strthinkick(r6,elem->mpole.PolynomA,elem->mpole.PolynomB,%STEP%,elem->mpole.MaxOrder)");

  // Pure quad
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("quadthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,%STEP%)");
  // Pure sextu
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("sextuthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,%STEP%)");
  // Pure octu
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("octuthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,%STEP%)");

  integrator.generateCode(code);

  generateQuadFringeExit(code,info);
  generateApertures(code,info);
  generateExit(code,info);

}

void StrMPoleSymplectic4Pass::generateQuadFringeEnter(std::string& code, PassMethodInfo *info) noexcept {

  if(info->doQuadEnter)
    code.append(
          "  if (elem->mpole.FringeQuadEntrance) {\n"
          "    quad_fringe(r6,elem->mpole.PolynomB[1],(AT_FLOAT)1.0,p_norm);\n"
          "  }\n"
  );

}

void StrMPoleSymplectic4Pass::generateQuadFringeExit(std::string& code, PassMethodInfo *info) noexcept {

  if(info->doQuadExit)
    code.append(
          "  if(elem->mpole.FringeQuadExit) {\n"
          "    quad_fringe(r6,elem->mpole.PolynomB[1],(AT_FLOAT)(-1.0),p_norm);\n"
          "  }\n"
  );

}

void StrMPoleSymplectic4Pass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  AbstractGPU::getInstance()->getDeviceFunctionQualifier(ftype);

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
          "  r[3] += L * (K * ((AT_FLOAT)2.0* r[0] * r[2]) + A0);\n"
          "}\n"
  );

  // Octu
  code.append(
          ftype +
          "void octuthinkick(AT_FLOAT* r,AT_FLOAT A0,AT_FLOAT B0,AT_FLOAT K,AT_FLOAT L) {\n"
          "  AT_FLOAT x2 = r[0]*r[0];\n"
          "  AT_FLOAT y2 = r[2]*r[2];\n"
          "  r[1] -= L * ((K * r[0] * (x2 - (AT_FLOAT)3.0*y2)) + B0);\n"
          "  r[3] += L * ((K * r[2] * ((AT_FLOAT)3.0*x2 - y2)) + A0);\n"
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
          "  AT_FLOAT u = p_norm * b2 / (AT_FLOAT)12.0;\n"
          "  AT_FLOAT x2 = r[0] * r[0];\n"
          "  AT_FLOAT z2 = r[2] * r[2];\n"
          "  AT_FLOAT xz = r[0] * r[2];\n"
          "  AT_FLOAT gx = u * (x2 + (AT_FLOAT)3.0 * z2) * r[0];\n"
          "  AT_FLOAT gz = u * (z2 + (AT_FLOAT)3.0 * x2) * r[2];\n"
          "  AT_FLOAT r1tmp = (AT_FLOAT)3.0 * u * ((AT_FLOAT)2.0 * xz * r[3] - (x2 + z2) * r[1]);\n"
          "  AT_FLOAT r3tmp = (AT_FLOAT)3.0 * u * ((AT_FLOAT)2.0 * xz * r[1] - (x2 + z2) * r[3]);\n"
          "  r[0] += sign*gx;\n"
          "  r[2] -= sign*gz;\n"
          "  r[5] -= sign * (gz * r[3] - gx * r[1]) * p_norm;\n"
          "  r[1] += sign*r1tmp;\n"
          "  r[3] -= sign*r3tmp;\n"
          "}\n"
  );

}
