#include "DriftPass.h"

DriftPass::DriftPass() noexcept : IdentityPass() {

}

void DriftPass::getParameters(AbstractInterface *param, PASSMETHOD_INFO *info) {

  // Retrieve param from super class
  IdentityPass::getParameters(param,info);
  Length = param->getOptionalDouble("Length",0);

}

void DriftPass::generateGPUKernel(std::string& code, PASSMETHOD_INFO *info) noexcept {

  code.append("__device__ void DriftPass(AT_FLOAT* r6,ELEMENT* elem) {\n");
  generateEnter(code,info);
  code.append("  AT_FLOAT p_norm = 1.0 / (1.0 + r6[4]);\n");
  code.append("  fastdrift(r6,elem->Length,p_norm);\n");
  generateExit(code,info);
  code.append("}\n");

}

void DriftPass::generateCall(std::string& code) noexcept {

  code.append("      case DRIFT:\n");
  code.append("        DriftPass(r6,elemPtr);\n");
  code.append("        break;\n");

}
