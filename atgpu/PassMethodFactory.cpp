#include "PassMethodFactory.h"
#include "IdentityPass.h"
#include "DriftPass.h"
#include <string.h>
using namespace std;

PassMethodFactory *PassMethodFactory::handler = nullptr;

void PassMethodFactory::reset() {

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

}

void PassMethodFactory::generatePassMethodsCalls(std::string& code) {
  code.append(callCode);
}

PassMethodFactory *PassMethodFactory::getInstance() {
  if( handler== nullptr )
    handler = new PassMethodFactory();
  return handler;
}

