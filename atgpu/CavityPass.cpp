#include "CavityPass.h"
#include "AbstractGPU.h"
#include "../atintegrators/atconstants.h"

using namespace std;

CavityPass::CavityPass() noexcept : IdentityPass() {
}

CavityPass::~CavityPass() noexcept {
}

void CavityPass::getParameters(AbstractInterface *param, PassMethodInfo *info) {

  // Retrieve param from super class
  IdentityPass::getParameters(param,info);

  elemData.Type = CAVITYPASS;
  elemData.Length = param->getDouble("Length");
  AT_FLOAT Energy = param->getDouble("Energy");
  AT_FLOAT Voltage = param->getDouble("Voltage");
  Frequency = param->getDouble("Frequency");
  elemData.NV = Voltage/Energy;
  elemData.FC = TWOPI * Frequency / C0;
  elemData.TimeLag = param->getOptionalDouble("TimeLag",0);
  elemData.PhaseLag = param->getOptionalDouble("PhaseLag",0);

}


void CavityPass::generateCode(std::string& code, PassMethodInfo *info,SymplecticIntegrator &integrator) noexcept {

  generateEnter(code,info);
  generateApertures(code,info);

  // Kick/Drift methods are defined in PassMethodFactory
  code.append(
          "  if(elem->Length == 0.0) {\n"
          "    r6[4] += -elem->NV*sin(elem->FC*(r6[5]-elem->TimeLag)-elem->PhaseLag);\n"
          "  } else {\n"
          "    fastdrift(r6,elem->Length/2.0,1.0/(1.0+r6[4]));\n"
          "    r6[4] += -elem->NV*sin(elem->FC*(r6[5]-elem->TimeLag)-elem->PhaseLag);\n"
          "    fastdrift(r6,elem->Length/2.0,1.0/(1.0+r6[4]));\n"
          "  }\n"
  );

  generateApertures(code,info);
  generateExit(code,info);

}

void CavityPass::generateUtilsFunction(std::string& code, PassMethodInfo *info) noexcept {

}
