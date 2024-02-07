#include "ExactMultipoleRadPass.h"
#include "../atintegrators/atconstants.h"
#include "AbstractGPU.h"
#include "PassMethodFactory.h"

using namespace std;

ExactMultipoleRadPass::ExactMultipoleRadPass() noexcept : ExactMultipolePass() {
}

ExactMultipoleRadPass::~ExactMultipoleRadPass() noexcept {
}

void ExactMultipoleRadPass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  ExactMultipolePass::getParameters(param,info);

  elemData.Type = EXACTMULTIPOLERADPASS;
  AT_FLOAT E0 = param->getDouble("Energy");
  elemData.mpole.CRAD =  CGAMMA*E0*E0*E0/(TWOPI*1e27);

}

void ExactMultipoleRadPass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateEnter(code,info);
  generateApertures(code,info);
  generateQuadFringeEnter(code,info);

  integrator.resetMethods();
  // Default straight element
  integrator.addDriftMethod("exact_drift(r6,%STEP%)");
  integrator.addKickMethod("exact_strthinkickrad(r6,elem->mpole.PolynomA,elem->mpole.PolynomB,%STEP%,elem->mpole.MaxOrder,"
                           "elem->mpole.CRAD)");
  // Pure quad
  integrator.addDriftMethod("exact_drift(r6,%STEP%)");
  integrator.addKickMethod("exact_quadthinkickrad(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,%STEP%,"
                           "elem->mpole.CRAD)");

  integrator.generateCode(code);

  // Convert absolute path length to path lengthening
  code.append("  r6[5] -= elem->Length;\n");

  generateQuadFringeExit(code,info);
  generateApertures(code,info);
  generateExit(code,info);

}

void ExactMultipoleRadPass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  AbstractGPU::getInstance()->getDeviceFunctionQualifier(ftype);


  string radCode;
  radCode.append(
          // 'denormalize' x and y prime
          "  AT_FLOAT p_norm = PNORM(r[4]);\n"
          "  AT_FLOAT xpr = r[1]*p_norm;\n"
          "  AT_FLOAT ypr = r[3]*p_norm;\n"
          // Compute normalized transverse field ( |B x v|^2 ) in straight element
          "  AT_FLOAT bx = ImSum;\n"
          "  AT_FLOAT by = ReSum;\n"
          "  AT_FLOAT B2P = SQR(bx) + SQR(by) - SQR(bx*xpr + by*ypr);\n"
          // Energy loss
          "  r[4] = r[4] - CRAD * SQR((AT_FLOAT)1+r[4]) * B2P * L * rsqrt((AT_FLOAT)1 - SQR(xpr) - SQR(ypr));\n"
          // Recalculate momentum from angles after losing energy for radiation
          "  r[1] = xpr + xpr * r[4];\n"
          "  r[3] = ypr + ypr * r[4];\n"
          // Kick
          "  r[1] -=  L*ReSum;\n"
          "  r[3] +=  L*ImSum;\n"
  );

  // Generic kick in straight element (rad)
  code.append(
          ftype +
          "void exact_strthinkickrad(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order,AT_FLOAT CRAD) {\n"
          + PassMethodFactory::polyLoop
          + radCode +
          "}\n"
  );

  // Generic kick in quad (rad)
  code.append(
          ftype +
          "void exact_quadthinkickrad(AT_FLOAT* r,AT_FLOAT A0,AT_FLOAT B0,AT_FLOAT K,AT_FLOAT L,AT_FLOAT CRAD) {\n"
          "  AT_FLOAT ReSum = (K * r[0] + B0);\n"
          "  AT_FLOAT ImSum = (K * r[2] + A0);\n"
          + radCode +
          "}\n"
  );

}