#include "AbstractInterface.h"

AbstractInterface *AbstractInterface::handler = nullptr;

void AbstractInterface::setHandler(AbstractInterface *obj) {
  handler = obj;
}

AbstractInterface *AbstractInterface::getInstance() {
  if( handler== nullptr )
    throw ("AbstractInterface: handler not set");
  return handler;
}
