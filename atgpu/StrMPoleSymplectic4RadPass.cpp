#include "StrMPoleSymplectic4RadPass.h"
#include "../atintegrators/atconstants.h"
#include "AbstractGPU.h"
#include <math.h>
#include "PassMethodFactory.h"

using namespace std;

StrMPoleSymplectic4RadPass::StrMPoleSymplectic4RadPass() noexcept : StrMPoleSymplectic4Pass() {

}

StrMPoleSymplectic4RadPass::~StrMPoleSymplectic4RadPass() noexcept {
}

void StrMPoleSymplectic4RadPass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  StrMPoleSymplectic4Pass::getParameters(param,info);

  if( elemData.Type==DRIFTPASS )
    // Do not change straight element retrograded to drift
    return;

  elemData.Type = STRMPOLESYMPLECTIC4RADPASS;
  AT_FLOAT E0 = param->getDouble("Energy");
  elemData.CRAD =  CGAMMA*E0*E0*E0/(TWOPI*1e27);

}

void StrMPoleSymplectic4RadPass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  code.append("  AT_FLOAT p_norm = 1.0 / (1.0 + r6[4]);\n");
  generateEnter(code,info);
  generateApertures(code,info);
  generateQuadFringeEnter(code,info);

  integrator.resetMethods();
  // Default straight element
  integrator.addDriftMethod("p_norm=1.0/(1.0 + r6[4]);fastdrift(r6,%STEP%,p_norm)");
  integrator.addKickMethod("strthinkickrad(r6,elem->PolynomA,elem->PolynomB,%STEP%,elem->MaxOrder,elem->CRAD,p_norm)");

  integrator.generateCode(code);

  if(integrator.getLastKickWeight()!=0.0)
    code.append("  p_norm = 1.0 / (1.0 + r6[4]);\n");

  generateQuadFringeExit(code,info);
  generateApertures(code,info);
  generateExit(code,info);

}

void StrMPoleSymplectic4RadPass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  getGPUFunctionQualifier(ftype);

  // Generic kick in straight element (rad)
  code.append(
          ftype +
          "void strthinkickrad(AT_FLOAT* r,const AT_FLOAT* A,const AT_FLOAT* B,AT_FLOAT L,int max_order,AT_FLOAT CRAD,AT_FLOAT p_norm) {\n"
          + PassMethodFactory::polyLoop +
          "  AT_FLOAT xpr = r[1]*p_norm;\n"
          "  AT_FLOAT ypr = r[3]*p_norm;\n"
          "  AT_FLOAT xpr2 = SQR(xpr);\n"
          "  AT_FLOAT ypr2 = SQR(ypr);\n"
          "  AT_FLOAT v_norm2 = 1.0/(1.0 + xpr2 + ypr2);\n"
          // Compute normalized transverse field ( |B x v|^2 ) in straight element
          "  AT_FLOAT bx = ImSum;\n"
          "  AT_FLOAT by = ReSum;\n"
          "  AT_FLOAT B2P = ( SQR(by) + SQR(bx) + SQR(bx*ypr - by*xpr) )*v_norm2;\n"
          // Energy loss
          "  r[4] = r[4] - CRAD*SQR(1.0+r[4])*B2P*(1.0 + (xpr2+ypr2)*0.5) * L;\n"
          // recalculate momentums from angles after losing energy for radiation
          "  r[1] = xpr * (1.0+r[4]);\n"
          "  r[3] = ypr * (1.0+r[4]);\n"
          // Kick
          "  r[1] -=  L*ReSum;\n"
          "  r[3] +=  L*ImSum;\n"
          "}\n"
  );

}