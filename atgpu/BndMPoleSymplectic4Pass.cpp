#include "BndMPoleSymplectic4Pass.h"
#include "AbstractGPU.h"
#include <math.h>
#include <iostream>

using namespace std;

BndMPoleSymplectic4Pass::BndMPoleSymplectic4Pass(SymplecticIntegrator& integrator) noexcept : StrMPoleSymplectic4Pass(integrator) {

}

BndMPoleSymplectic4Pass::~BndMPoleSymplectic4Pass() noexcept {
}

void BndMPoleSymplectic4Pass::getParameters(AbstractInterface *param, PASSMETHOD_INFO *info) {

  AT_FLOAT sedge;
  AT_FLOAT cedge;
  AT_FLOAT tedge;

  // Retrieve param from super class
  StrMPoleSymplectic4Pass::getParameters(param,info);

  elemData.Type = BEND;
  elemData.irho = param->getDouble("BendingAngle") / elemData.Length;
  AT_FLOAT gap = param->getOptionalDouble("FullGap", 0);

  elemData.FringeBendEntrance = param->getOptionalInt("FringeBendEntrance",1);
  AT_FLOAT fint = param->getOptionalDouble("FringeInt1", 0);
  elemData.EntranceAngle = param->getOptionalDouble("EntranceAngle", 0);
  sedge = sin(elemData.EntranceAngle);
  cedge = cos(elemData.EntranceAngle);
  tedge = tan(elemData.EntranceAngle);
  elemData.FringeCorrEntranceX = elemData.irho*tedge;
  elemData.FringeCorrEntranceY = elemData.irho*gap*fint*(1.0+sedge*sedge)/cedge;

  elemData.FringeBendExit = param->getOptionalInt("FringeBendExit",1);
  fint = param->getOptionalDouble("FringeInt2", 0);
  elemData.ExitAngle = param->getOptionalDouble("ExitAngle", 0);
  sedge = sin(elemData.ExitAngle);
  cedge = cos(elemData.ExitAngle);
  tedge = tan(elemData.ExitAngle);
  elemData.FringeCorrExitX = elemData.irho*tedge;
  elemData.FringeCorrExitY = elemData.irho*gap*fint*(1.0+sedge*sedge)/cedge;

}


void BndMPoleSymplectic4Pass::generateGPUKernel(std::string& code, PASSMETHOD_INFO *info,SymplecticIntegrator& integrator) noexcept {

  AbstractGPU *gpu = AbstractGPU::getInstance();
  string ftype;
  gpu->getDeviceFunctionQualifier(ftype);
  if(!ftype.empty()) ftype.append(" ");

  code.append( ftype + "void BndMPoleSymplectic4Pass(AT_FLOAT* r6,ELEMENT* elem) {\n");

  code.append(
          "  AT_FLOAT p_norm = 1.0 / (1.0 + r6[4]);\n"
  );

  generateEnter(code,info);
  generateApertures(code,info);
  generateKickAngle(code,info);
  generateBendFringeEnter(code,info);
  generateQuadFringeEnter(code,info);
  integrator.generateCode(code,"elem->NumIntSteps","fastdrift","bndthinkick","elem->irho");
  generateQuadFringeExit(code,info);
  generateBendFringeExit(code,info);
  generateKickAngleRestore(code,info);
  generateApertures(code,info);
  generateExit(code,info);
  code.append("}\n");

}

void BndMPoleSymplectic4Pass::generateCall(std::string& code) noexcept {

  code.append(
          "      case BEND:\n"
          "        BndMPoleSymplectic4Pass(r6,elemPtr);\n"
          "        break;\n"
  );

}

void BndMPoleSymplectic4Pass::generateBendFringeEnter(std::string& code, PASSMETHOD_INFO *info) noexcept {
  code.append("  edge_fringe(r6,p_norm,elem->EntranceAngle,elem->irho,elem->FringeCorrEntranceX,elem->FringeCorrEntranceY,elem->FringeBendEntrance);\n");
}

void BndMPoleSymplectic4Pass::generateBendFringeExit(std::string& code, PASSMETHOD_INFO *info) noexcept {
  code.append("  edge_fringe(r6,p_norm,elem->ExitAngle,elem->irho,elem->FringeCorrExitX,elem->FringeCorrExitY,elem->FringeBendExit);\n");
}


