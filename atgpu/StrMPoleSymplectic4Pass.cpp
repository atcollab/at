#include "StrMPoleSymplectic4Pass.h"
#include <vector>
#include <string.h>
#include <math.h>

using namespace std;

StrMPoleSymplectic4Pass::StrMPoleSymplectic4Pass() noexcept : IdentityPass() {

}

void StrMPoleSymplectic4Pass::getParameters(AbstractInterface *param, PASSMETHOD_INFO *info) {

  // Retrieve param from super class
  IdentityPass::getParameters(param,info);

  elemData.Type = MPOLE;
  elemData.Length = param->getDouble("Length");
  elemData.MaxOrder = param->getInt("MaxOrder");
  getPolynom(elemData.PolynomA,param,"PolynomA",elemData.MaxOrder);
  getPolynom(elemData.PolynomB,param,"PolynomB",elemData.MaxOrder);
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
  elemData.FringeQuadEntrance = param->getOptionalDouble("FringeQuadEntrance", 0);
  elemData.FringeQuadExit = param->getOptionalDouble("FringeQuadExit", 0);
  if( elemData.MaxOrder>=1 && elemData.PolynomB[1]==0.0 ) {
    // No quad strength
    elemData.FringeQuadEntrance = false;
    elemData.FringeQuadExit = false;
  }

  info->doQuadEnter |= elemData.FringeQuadEntrance;
  info->doQuadExit |= elemData.FringeQuadExit;

}

void StrMPoleSymplectic4Pass::getPolynom(AT_FLOAT *dest,AbstractInterface *param,const std::string& name,int order) {

  if( order==0 )
    return;

  static vector<int64_t> shape;
  AT_FLOAT *P = param->getNativeDoubleArray(name,shape);
  if(shape.size()!=1 || shape[0]<elemData.MaxOrder)
    throw string(name + ", wrong dimension: (" + to_string(elemData.MaxOrder) + ") expected gut got " +
                 AbstractInterface::getShapeStr(shape));
  memcpy(dest, P, elemData.MaxOrder*sizeof(AT_FLOAT));

}

void StrMPoleSymplectic4Pass::generateGPUKernel(std::string& code, PASSMETHOD_INFO *info) noexcept {

  code.append("__device__ void MPoleSymplectic4Pass(AT_FLOAT* r6,ELEMENT* elem) {\n");
  generateEnter(code,info);
  code.append("  AT_FLOAT p_norm = 1.0 / (1.0 + r6[4]);\n");
  code.append("  fastdrift(r6,elem->Length,p_norm);\n");
  generateExit(code,info);
  code.append("}\n");

}

void StrMPoleSymplectic4Pass::generateCall(std::string& code) noexcept {

  code.append("      case DRIFT:\n");
  code.append("        MPoleSymplectic4Pass(r6,elemPtr);\n");
  code.append("        break;\n");

}
