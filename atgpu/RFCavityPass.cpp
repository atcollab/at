#include "RFCavityPass.h"
#include "AbstractGPU.h"
#include "../atintegrators/atconstants.h"

using namespace std;

RFCavityPass::RFCavityPass() noexcept : CavityPass() {
}

RFCavityPass::~RFCavityPass() noexcept {
}

void RFCavityPass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  CavityPass::getParameters(param,info);

  elemData.Type = RFCAVITYPASS;
  elemData.Length = param->getDouble("Length");

}

void RFCavityPass::postInit(RING_PARAM *param) {

  // Harmonic number correction
  double T0 = param->Length / C0;
  elemData.HC = C0*( round( Frequency * T0 ) / Frequency  - T0 );

}

void RFCavityPass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateEnter(code,info);
  generateApertures(code,info);

  // Kick/Drift methods are defined in PassMethodFactory
  code.append(
          "  if(elem->Length == 0.0) {\n"
          "    r6[4] += -elem->NV*sin(elem->FC*(r6[5] - elem->TimeLag - elem->HC*turn) - elem->PhaseLag);\n"
          "  } else {\n"
          "    fastdrift(r6,elem->Length/2.0,1.0/(1.0+r6[4]));\n"
          "    r6[4] += -elem->NV*sin(elem->FC*(r6[5] - elem->TimeLag - elem->HC*turn) - elem->PhaseLag);\n"
          "    fastdrift(r6,elem->Length/2.0,1.0/(1.0+r6[4]));\n"
          "  }\n"
  );

  generateApertures(code,info);
  generateExit(code,info);

}

void RFCavityPass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

}
