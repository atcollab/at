#include "ExactMultipolePass.h"
#include "PassMethodFactory.h"
#include "AbstractGPU.h"

using namespace std;

ExactMultipolePass::ExactMultipolePass() noexcept : StrMPoleSymplectic4Pass() {
}

ExactMultipolePass::~ExactMultipolePass() noexcept {
}

void ExactMultipolePass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  StrMPoleSymplectic4Pass::getParameters(param,info);

  elemData.Type = EXACTMULTIPOLEPASS;

}


void ExactMultipolePass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateEnter(code,info);
  generateApertures(code,info);
  generateQuadFringeEnter(code,info);

  integrator.resetMethods();
  // Default straight magnet
  integrator.addDriftMethod("exact_drift(r6,%STEP%)");
  integrator.addKickMethod("strthinkick(r6,elem->mpole.PolynomA,elem->mpole.PolynomB,%STEP%,elem->mpole.MaxOrder)");

  // Pure quad
  integrator.addDriftMethod("exact_drift(r6,%STEP%)");
  integrator.addKickMethod("quadthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,%STEP%)");
  // Pure sextu
  integrator.addDriftMethod("exact_drift(r6,%STEP%)");
  integrator.addKickMethod("sextuthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,%STEP%)");
  // Pure octu
  integrator.addDriftMethod("exact_drift(r6,%STEP%)");
  integrator.addKickMethod("octuthinkick(r6,elem->mpole.PolynomA[0],elem->mpole.PolynomB[0],elem->mpole.K,%STEP%)");

  integrator.generateCode(code);

  // Convert absolute path length to path lengthening
  code.append("  r6[5] -= elem->Length;\n");

  generateQuadFringeExit(code,info);
  generateApertures(code,info);
  generateExit(code,info);

}

void ExactMultipolePass::generateQuadFringeEnter(std::string& code, PassMethodInfo *info) noexcept {

  if(info->doQuadEnter)
    code.append(
            "  if (elem->mpole.FringeQuadEntrance) {\n"
            "    multipole_fringe(r6,elem->mpole.PolynomA,elem->mpole.PolynomB,elem->mpole.MaxOrder,(AT_FLOAT)1.0,false);\n"
            "  }\n"
    );

}

void ExactMultipolePass::generateQuadFringeExit(std::string& code, PassMethodInfo *info) noexcept {

  if(info->doQuadExit)
    code.append(
            "  if(elem->mpole.FringeQuadExit) {\n"
            "    multipole_fringe(r6,elem->mpole.PolynomA,elem->mpole.PolynomB,elem->mpole.MaxOrder,(AT_FLOAT)-1.0,false);\n"
            "  }\n"
    );

}



void ExactMultipolePass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

  string ftype;
  AbstractGPU::getInstance()->getDeviceFunctionQualifier(ftype);

  // PTC multipole_fringer
  // Forest 13.29
  // not re-derived and checked
  // note this is the sum over n of Forest 13.29
  // one for each multipole component
  code.append(
          ftype +
          "void multipole_fringe(AT_FLOAT *r, AT_FLOAT *A, AT_FLOAT *B, int max_order,AT_FLOAT sign, int skip_b0) {\n"
          "  AT_FLOAT FX = (AT_FLOAT)0;\n"
          "  AT_FLOAT FY = (AT_FLOAT)0;\n"
          "  AT_FLOAT FX_X = (AT_FLOAT)0;\n"
          "  AT_FLOAT FX_Y = (AT_FLOAT)0;\n"
          "  AT_FLOAT FY_X = (AT_FLOAT)0;\n"
          "  AT_FLOAT FY_Y = (AT_FLOAT)0;\n"
          "  AT_FLOAT RX = (AT_FLOAT)1;\n"
          "  AT_FLOAT IX = (AT_FLOAT)0;\n"
          // Invariant is: (j is the index, i is the complex unit)
          // RX+IXi = (x + iy)^j
          "  for (int n = 0; n <= max_order; n++) {\n"
          "    AT_FLOAT j = (AT_FLOAT)(n + 1);\n"
          //   Complex mult
          "    AT_FLOAT DRX = RX;\n"
          "    AT_FLOAT DIX = IX;\n"
          "    RX = DRX * r[0] - DIX * r[2];\n"
          "    IX = DRX * r[2] + DIX * r[0];\n"
          "    AT_FLOAT f1 = -sign * (AT_FLOAT)0.25 / (j + (AT_FLOAT)1);\n"
          "    AT_FLOAT nf = (j + (AT_FLOAT)2) / j;\n"
          "    AT_FLOAT Bn = (n == 0 && skip_b0)?(AT_FLOAT)0:B[n];\n"
          "    AT_FLOAT An = A[n];\n"
          "    AT_FLOAT U = f1 * (Bn * RX - An * IX);\n"
          "    AT_FLOAT V = f1 * (Bn * IX + An * RX);\n"
          "    AT_FLOAT DU = j * f1 * (Bn * DRX - An * DIX);\n"
          "    AT_FLOAT DV = j * f1 * (Bn * DIX + An * DRX);\n"
          "    FX += U * r[0] + nf * V * r[2];\n"
          "    FY += U * r[2] - nf * V * r[0];\n"
          "    FX_X +=  DU * r[0] + U      + nf * r[2] * DV;\n"
          "    FY_X +=  DU * r[2] - V * nf - nf * r[0] * DV;\n"
          "    FX_Y += -DV * r[0] + V * nf + nf * r[2] * DU;\n"
          "    FY_Y += -DV * r[2] + U      - nf * r[0] * DU;\n"
          "  }\n"
          "  AT_FLOAT DEL = PNORM(r[4]);\n"
          // Solve 2x2 matrix equation (Cramer's rule) for px and py
          "  AT_FLOAT M00 = (AT_FLOAT)1 - FX_X * DEL;\n"
          "  AT_FLOAT M10 = -FY_X * DEL;\n"
          "  AT_FLOAT M01 = -FX_Y * DEL;\n"
          "  AT_FLOAT M11 = (AT_FLOAT)1 - FY_Y * DEL;\n"
          "  AT_FLOAT idet = (AT_FLOAT)1 / (M00*M11-M10*M01);\n"
          "  AT_FLOAT px = (M11 * r[1] - M10 * r[3]) * idet;\n"
          "  AT_FLOAT py = (M00 * r[3] - M01 * r[1]) * idet;\n"
          "  r[0] = r[0] - FX * DEL;\n"
          "  r[1] = px;\n"
          "  r[2] = r[2] - FY * DEL;\n"
          "  r[3] = py;\n"
          "  r[5] = r[5] - (r[1] * FX + r[3] * FY) * DEL * DEL;\n"
          "}\n"
  );

}
