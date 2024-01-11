#include "AbstractInterface.h"
#include <string>
#include <string.h>

using namespace std;


AbstractInterface *AbstractInterface::handler = nullptr;

void AbstractInterface::setHandler(AbstractInterface *obj) {
  handler = obj;
}

AbstractInterface *AbstractInterface::getInstance() {
  if( handler== nullptr )
    throw string("AbstractInterface: handler not set");
  return handler;
}

AT_FLOAT *AbstractInterface::getDoubleArray(const std::string& name, std::vector<int64_t> expectedShape) {

  std::vector<int64_t> shape;
  AT_FLOAT *array = getNativeDoubleArray(name,shape);

  if(shape.size()!=expectedShape.size()) {
    throw string(name + " array attribute has shape "+ getShapeStr(shape) +" but " +
                 getShapeStr(expectedShape) + " expected");
  }

  bool ok=true;
  size_t d = 0;
  while(ok && d<expectedShape.size()) {
    ok = shape[d] == expectedShape[d];
    d++;
  }
  if( !ok ) {
    throw string(name + " array attribute has shape "+ getShapeStr(shape) +" but " +
                 getShapeStr(expectedShape) + " expected");
  }

  return array;

}

void AbstractInterface::getOptional1DArray(AT_FLOAT **dest,const std::string& name,int length) {

  try {
    get1DArray(dest,name,length);
  } catch (string&) {
  }

}

void AbstractInterface::get1DArray(AT_FLOAT **dest,const std::string& name,int length) {

  *dest = nullptr;

  static vector<int64_t> shape;
  AT_FLOAT *P = getNativeDoubleArray(name,shape);
  if(shape.size()!=1 || shape[0]<length) {
    throw string(name + ", wrong dimension: (" + to_string(length) + ") expected gut got " +
                 AbstractInterface::getShapeStr(shape));
  }

  *dest = P;

}

std::string AbstractInterface::getShapeStr(std::vector<int64_t>& shape) {
  std::string ret = "(";
  for(size_t i=0;i<shape.size();i++) {
    ret.append(std::to_string(shape[i]));
    if(i<shape.size()-1) ret.append("x");
  }
  ret.append(")");
  return ret;
}

void AbstractInterface::getShapeFromStr(std::vector<int64_t>& shape,std::string& str) {

  if(str.size()<3 || str[0]!='(' || str[str.size()-1]!=')')
    throw string("Invalid shape definition, got " + str + " expected (N[xMx...])");

  vector<string> stokens;
  split(stokens, str, 'x');
  stokens[0].erase(stokens[0].begin());
  stokens[stokens.size()-1].pop_back();
  for (const auto &stoken: stokens)
    shape.push_back(stoi(stoken));

}

void AbstractInterface::split(vector<string> &tokens, const string &text, char sep) {

  size_t start = 0, end = 0;
  tokens.clear();

  while ((end = text.find(sep, start)) != string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }

  tokens.push_back(text.substr(start));

}

int AbstractInterface::getOptionalInt(const std::string& name,int defaultValue) {

  try {
    return getInt(name);
  } catch (string&) {
    return defaultValue;
  }

}

AT_FLOAT AbstractInterface::getOptionalDouble(const std::string& name, AT_FLOAT defaultValue) {

  try {
    return getDouble(name);
  } catch (string&) {
    return defaultValue;
  }

}

AT_FLOAT *AbstractInterface::getOptionalDoubleArray(const std::string &name, std::vector<int64_t> expectedShape) {

  try {
    return getDoubleArray(name,expectedShape);
  } catch (string&) {
    return nullptr;
  }

}
