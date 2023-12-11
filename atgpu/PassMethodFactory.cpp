#include "PassMethodFactory.h"
#include "IdentityPass.h"
#include "DriftPass.h"
#include "StrMPoleSymplectic4Pass.h"
#include <string.h>

using namespace std;

PassMethodFactory::PassMethodFactory(SymplecticIntegrator& integrator) noexcept:integrator(integrator) {

  // Flags for pass methods
  memset(passMethodInfos, 0, sizeof(passMethodInfos));

}

AbstractElement *PassMethodFactory::createElement(std::string& passMethod) {

  AbstractElement *elem;

  // Get a pointer to the abstract interface
  AbstractInterface *I = AbstractInterface::getInstance();

  if(passMethod=="IdentityPass") {
    elem = new IdentityPass();
    elem->getParameters(I,&passMethodInfos[IDENTITY]);
  } else if (passMethod=="DriftPass") {
    elem = new DriftPass();
    elem->getParameters(I,&passMethodInfos[DRIFT]);
  } else if (passMethod=="StrMPoleSymplectic4Pass") {
    elem = new StrMPoleSymplectic4Pass(integrator);
    elem->getParameters(I,&passMethodInfos[MPOLE]);
  } else {
    throw string("Not implemented PassMethod: " + passMethod);
  }

  return elem;

}

void PassMethodFactory::generatePassMethods(std::string& code) {

  callCode.clear();

  if( passMethodInfos[IDENTITY].used ) {
    IdentityPass::generateGPUKernel(code,&passMethodInfos[IDENTITY]);
    IdentityPass::generateCall(callCode);
  }
  if( passMethodInfos[DRIFT].used ) {
    DriftPass::generateGPUKernel(code,&passMethodInfos[DRIFT]);
    DriftPass::generateCall(callCode);
  }
  if( passMethodInfos[MPOLE].used ) {
    StrMPoleSymplectic4Pass::generateGPUKernel(code,&passMethodInfos[MPOLE],integrator);
    StrMPoleSymplectic4Pass::generateCall(callCode);
  }

}

void PassMethodFactory::generatePassMethodsCalls(std::string& code) {
  code.append(callCode);
}

