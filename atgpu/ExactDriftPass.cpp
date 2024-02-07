#include "ExactDriftPass.h"

using namespace std;

ExactDriftPass::ExactDriftPass() noexcept : DriftPass() {}

void ExactDriftPass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  DriftPass::getParameters(param,info);

  elemData.Type = EXACTDRIFTPASS;

}

void ExactDriftPass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateEnter(code,info);
  generateApertures(code,info);
  code.append("  exact_drift(r6,elem->Length);\n");
  // Convert absolute path length to path lengthening
  code.append("  r6[5] -= elem->Length;\n");
  generateApertures(code,info);
  generateExit(code,info);

}

void ExactDriftPass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  AbstractGPU::getInstance()->getDeviceFunctionQualifier(ftype);

  //Exact drift
  code.append(
          ftype +
          "void exact_drift(AT_FLOAT* r,AT_FLOAT L) {\n"
          "  AT_FLOAT p_norm = L * rsqrt(SQR(1 + r[4]) - SQR(r[1]) - SQR(r[3]));\n"
          "  r[0] += p_norm * r[1];\n"
          "  r[2] += p_norm * r[3];\n"
          // Absolute path length
          "  r[5] += p_norm * (1.0 + r[4]);\n"
          "}\n"
  );

}
