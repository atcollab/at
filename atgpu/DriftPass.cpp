#include "DriftPass.h"

using namespace std;

DriftPass::DriftPass() noexcept : IdentityPass() {}

void DriftPass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  IdentityPass::getParameters(param,info);

  elemData.Type = DRIFTPASS;
  elemData.Length = param->getDouble("Length");

}

void DriftPass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateEnter(code,info);
  generateApertures(code,info);
  code.append("  AT_FLOAT p_norm = PNORM(r6[4]);\n");
  code.append("  fastdrift(r6,elem->Length,p_norm);\n");
  generateApertures(code,info);
  generateExit(code,info);

}

void DriftPass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  AbstractGPU::getInstance()->getDeviceFunctionQualifier(ftype);

  //Drift (small angle)
  code.append(
          ftype +
          "void fastdrift(AT_FLOAT* r,AT_FLOAT L,AT_FLOAT p_norm) {\n"
          "  r[0] += p_norm * L * r[1];\n"
          "  r[2] += p_norm * L * r[3];\n"
          "  r[5] += p_norm * p_norm * L * (r[1] * r[1] + r[3] * r[3]) * (AT_FLOAT)0.5;\n"
          "}\n"
  );

}
