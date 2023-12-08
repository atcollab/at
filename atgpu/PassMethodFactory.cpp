#include "PassMethodFactory.h"
#include "IdentityPass.h"
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
  } else {
    throw string("not implemented PassMethod: \" + passMethod");
  }

  return elem;

}

void PassMethodFactory::generatePassMethods(std::string& code) {

  callCode.clear();

  for(auto & i : passMethodInfos) {
    if( i.used ) {
      // IdentyPass is the super class
      IdentityPass::generateGPUKernel(code,&i);
      IdentityPass::generateCall(callCode);
    }
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

