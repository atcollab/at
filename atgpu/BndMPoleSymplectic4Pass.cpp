#include "BndMPoleSymplectic4Pass.h"
#include "AbstractGPU.h"
#include "PassMethodFactory.h"
#include <math.h>

using namespace std;

BndMPoleSymplectic4Pass::BndMPoleSymplectic4Pass() noexcept : StrMPoleSymplectic4Pass() {
}

BndMPoleSymplectic4Pass::~BndMPoleSymplectic4Pass() noexcept {
}

void BndMPoleSymplectic4Pass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  AT_FLOAT sedge;
  AT_FLOAT cedge;
  AT_FLOAT tedge;

  // Retrieve param from super class
  StrMPoleSymplectic4Pass::getParameters(param,info);

  elemData.Type = BNDMPOLESYMPLECTIC4PASS;
  elemData.mpole.irho = param->getDouble("BendingAngle") / elemData.Length;
  AT_FLOAT gap = param->getOptionalDouble("FullGap", 0);

  elemData.mpole.FringeBendEntrance = param->getOptionalInt("FringeBendEntrance",1);
  AT_FLOAT fint = param->getOptionalDouble("FringeInt1", 0);
  elemData.mpole.EntranceAngle = param->getOptionalDouble("EntranceAngle", 0);
  sedge = sin(elemData.mpole.EntranceAngle);
  cedge = cos(elemData.mpole.EntranceAngle);
  tedge = tan(elemData.mpole.EntranceAngle);
  elemData.mpole.FringeCorrEntranceX = elemData.mpole.irho*tedge;
  elemData.mpole.FringeCorrEntranceY = elemData.mpole.irho*gap*fint*(1.0+sedge*sedge)/cedge;

  elemData.mpole.FringeBendExit = param->getOptionalInt("FringeBendExit",1);
  fint = param->getOptionalDouble("FringeInt2", 0);
  elemData.mpole.ExitAngle = param->getOptionalDouble("ExitAngle", 0);
  sedge = sin(elemData.mpole.ExitAngle);
  cedge = cos(elemData.mpole.ExitAngle);
  tedge = tan(elemData.mpole.ExitAngle);
  elemData.mpole.FringeCorrExitX = elemData.mpole.irho*tedge;
  elemData.mpole.FringeCorrExitY = elemData.mpole.irho*gap*fint*(1.0+sedge*sedge)/cedge;

}


void BndMPoleSymplectic4Pass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  code.append("  AT_FLOAT p_norm = PNORM(r6[4]);\n");
  generateEnter(code,info);
  generateApertures(code,info);
  generateBendFringeEnter(code,info);
  generateQuadFringeEnter(code,info);

  integrator.resetMethods();
  // Default bend
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("bndthinkick(r6,elem->mpole.PolynomA,elem->mpole.PolynomB,%STEP%,elem->mpole.MaxOrder,elem->irho)");
  // Pure bend
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("bndthinkick0(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],%STEP%,elem->mpole.irho)");
  // Pure DQ
  integrator.addDriftMethod("fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("dqthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,%STEP%,elem->mpole.irho)");

  integrator.generateCode(code);

  generateQuadFringeExit(code,info);
  generateBendFringeExit(code,info);
  generateApertures(code,info);
  generateExit(code,info);

}

void BndMPoleSymplectic4Pass::fillGPUMemory(uint8_t *startAdd,ELEMENT *elemMem,uint64_t *offset) {

  StrMPoleSymplectic4Pass::fillGPUMemory(startAdd,elemMem,offset);

  // Update Subtype for bending
  if (isBending()) {
    elemData.SubType = 1;
  } else if (isDQ()) {
    elemData.SubType = 2;
  } else {
    elemData.SubType = 0;
  }
  elemMem->SubType = elemData.SubType;

}

void BndMPoleSymplectic4Pass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  AbstractGPU::getInstance()->getDeviceFunctionQualifier(ftype);

  // Bending fringe field correction (edge angle focusing)
  code.append(
          ftype +
          "void edge_fringe(AT_FLOAT* r,AT_FLOAT p_norm,AT_FLOAT edge_angle,AT_FLOAT irho,AT_FLOAT f_corrx,AT_FLOAT f_corry,uint32_t method) {\n"
          "  AT_FLOAT fy;\n"
          "  switch(method) {\n"
          "  case 0:\n"
          "    fy = f_corry;\n"
          "    break;\n"
          "  case 1:\n"
          "    fy = irho * tan(edge_angle - f_corry * p_norm);\n"
          "    break;\n"
          "  case 2:\n"
          "    fy = irho * tan(edge_angle - f_corry * p_norm) * p_norm;\n"
          "    break;\n"
          "  case 3:\n"
          "    fy = irho * tan(edge_angle - f_corry + r[1] * p_norm);\n"
          "    break;\n"
          "  }\n"
          "  r[1] += r[0] * f_corrx;\n"
          "  r[3] -= r[2] * fy;\n"
          "}\n"
  );

  // Generic kick in bending element
  code.append(
          ftype +
          "void bndthinkick(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order,AT_FLOAT irho) {\n"
          + PassMethodFactory::polyLoop +
          "  r[1] -= L * (ReSum - (r[4] - r[0] * irho) * irho);\n"
          "  r[3] += L * ImSum;\n"
          "  r[5] += L * irho * r[0];\n"
          "}\n"
  );

  // Pure bending
  code.append(
          ftype +
          "void bndthinkick0(AT_FLOAT* r,AT_FLOAT A0,AT_FLOAT B0,AT_FLOAT L,AT_FLOAT irho) {\n"
          "  r[1] -= L * (B0 - (r[4] - r[0] * irho) * irho);\n"
          "  r[3] += L * A0;\n"
          "  r[5] += L * irho * r[0];\n"
          "}\n"
  );

  // Pure DQ
  code.append(
          ftype +
          "void dqthinkick(AT_FLOAT* r,AT_FLOAT A0,AT_FLOAT B0,AT_FLOAT K,AT_FLOAT L,AT_FLOAT irho) {\n"
          "  r[1] -= L * ((K * r[0] + B0) - (r[4] - r[0] * irho) * irho);\n"
          "  r[3] += L * (K * r[2] + A0);\n"
          "  r[5] += L * irho * r[0];\n"
          "}\n"
  );

}

void BndMPoleSymplectic4Pass::generateBendFringeEnter(std::string& code, PassMethodInfo *info) noexcept {
  code.append("  edge_fringe(r6,p_norm,elem->mpole.EntranceAngle,elem->mpole.irho,elem->mpole.FringeCorrEntranceX,"
                  "              elem->mpole.FringeCorrEntranceY,elem->mpole.FringeBendEntrance);\n");
}

void BndMPoleSymplectic4Pass::generateBendFringeExit(std::string& code, PassMethodInfo *info) noexcept {
  code.append("  edge_fringe(r6,p_norm,elem->mpole.ExitAngle,elem->mpole.irho,elem->mpole.FringeCorrExitX,"
              "                  elem->mpole.FringeCorrExitY,elem->mpole.FringeBendExit);\n");
}

bool BndMPoleSymplectic4Pass::isBending() {
  return elemData.mpole.MaxOrder==0 ||
        (elemData.mpole.MaxOrder==1 && PolynomA[1]==0.0 && PolynomB[1]==0.0);
}

bool BndMPoleSymplectic4Pass::isDQ() {
  return elemData.mpole.MaxOrder==1 && PolynomA[1]==0.0;
}
