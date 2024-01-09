#include "BndMPoleSymplectic4RadPass.h"
#include "../atintegrators/atconstants.h"
#include "AbstractGPU.h"
#include "PassMethodFactory.h"
#include <math.h>

using namespace std;

BndMPoleSymplectic4RadPass::BndMPoleSymplectic4RadPass() noexcept : BndMPoleSymplectic4Pass() {
}

BndMPoleSymplectic4RadPass::~BndMPoleSymplectic4RadPass() noexcept {
}

void BndMPoleSymplectic4RadPass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  BndMPoleSymplectic4Pass::getParameters(param,info);

  elemData.Type = BNDMPOLESYMPLECTIC4RADPASS;
  AT_FLOAT E0 = param->getDouble("Energy");
  elemData.CRAD =  CGAMMA*E0*E0*E0/(TWOPI*1e27);

}

void BndMPoleSymplectic4RadPass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  code.append("  AT_FLOAT p_norm = 1.0 / (1.0 + r6[4]);\n");
  generateEnter(code,info);
  generateApertures(code,info);
  generateBendFringeEnter(code,info);
  generateQuadFringeEnter(code,info);

  integrator.resetMethods();
  // Default bend
  integrator.addDriftMethod("p_norm=1.0/(1.0 + r6[4]);fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("bndthinkickrad(r6,elem->PolynomA,elem->PolynomB,%STEP%,elem->MaxOrder,elem->irho,elem->CRAD,p_norm)");

  integrator.generateCode(code);

  if(integrator.getLastKickWeight()!=0.0)
    code.append("  p_norm = 1.0 / (1.0 + r6[4]);\n");

  generateQuadFringeExit(code,info);
  generateBendFringeExit(code,info);
  generateApertures(code,info);
  generateExit(code,info);

}

void BndMPoleSymplectic4RadPass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  getGPUFunctionQualifier(ftype);

  // Generic kick in bending element (rad)
  code.append(
          ftype +
          "void bndthinkickrad(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order,AT_FLOAT irho,AT_FLOAT CRAD,AT_FLOAT p_norm) {\n"
          + PassMethodFactory::polyLoop +
          "  AT_FLOAT xpr = r[1]*p_norm;\n"
          "  AT_FLOAT ypr = r[3]*p_norm;\n"
          "  AT_FLOAT xpr2 = SQR(xpr);\n"
          "  AT_FLOAT ypr2 = SQR(ypr);\n"
          "  AT_FLOAT oirho = 1.0 + r[0]*irho;\n"
          "  AT_FLOAT v_norm2 = 1.0/(SQR(oirho) + xpr2 + ypr2);\n"
          // Compute normalized transverse field ( |B x v|^2 ) in bending element
          "  AT_FLOAT bx = ImSum;\n"
          "  AT_FLOAT by = ReSum + irho;\n"
          "  AT_FLOAT B2P = ( SQR(by*oirho) + SQR(bx*oirho) + SQR(bx*ypr - by*xpr) )*v_norm2;\n"
          "  AT_FLOAT dp_0 = r[4];\n"
          // Energy loss
          "  r[4] = r[4] - CRAD*SQR(1.0+r[4])*B2P*(oirho + (xpr2+ypr2)*0.5) * L;\n"
          // recalculate momentums from angles after losing energy for radiation
          "  r[1] = xpr * (1.0+r[4]);\n"
          "  r[3] = ypr * (1.0+r[4]);\n"
          // Kick
          "  r[1] -=  L*(ReSum-(dp_0-r[0]*irho)*irho);\n"
          "  r[3] +=  L*ImSum;\n"
          // Path length
          "  r[5] +=  L*irho*r[0];\n"
          "}\n"
  );

}