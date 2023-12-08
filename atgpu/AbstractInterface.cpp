#include "AbstractInterface.h"
#include <string>

AbstractInterface *AbstractInterface::handler = nullptr;

void AbstractInterface::setHandler(AbstractInterface *obj) {
  handler = obj;
}

AbstractInterface *AbstractInterface::getInstance() {
  if( handler== nullptr )
    throw std::string("AbstractInterface: handler not set");
  return handler;
}

std::string AbstractInterface::getShapeStr(std::vector<int64_t>& shape) {
  std::string ret = "(";
  for(size_t i=0;i<shape.size();i++) {
    ret.append(std::to_string(shape[i]));
    if(i<shape.size()-1) ret.append(",");
  }
  ret.append(")");
  return ret;
}
