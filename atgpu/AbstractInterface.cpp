#include "AbstractInterface.h"
#include <string>

using namespace std;


AbstractInterface *AbstractInterface::handler = nullptr;

void AbstractInterface::setHandler(AbstractInterface *obj) {
  handler = obj;
}

bool AbstractInterface::isValidHandler() {
  return handler != nullptr;
}

AbstractInterface *AbstractInterface::getInstance() {
  if( handler== nullptr )
    throw string("AbstractInterface: handler not set");
  return handler;
}

void AbstractInterface::getOptional1DArray(double **dest,const std::string& name,int length) {

  try {
    get1DArray(dest,name,length);
  } catch (string&) {
  }

}

void AbstractInterface::getOptional1DArray(double **dest,const std::string& name,int *length) {

  *dest = nullptr;
  *length = 0;
  try {
    static vector<int64_t> shape;
    double *P = getNativeDoubleArray(name,shape);
    if(shape.size()!=1) {
      throw string(name + ", wrong dimension: signle dimension array expected gut got " +
                   AbstractInterface::getShapeStr(shape));
    }
    *dest = P;
    *length = (int)shape[0];
  } catch (string&) {
  }


}

void AbstractInterface::get1DArray(double **dest,const std::string& name,int length) {

  *dest = nullptr;

  static vector<int64_t> shape;
  double *P = getNativeDoubleArray(name,shape);
  if(shape.size()!=1 || shape[0]<length) {
    throw string(name + ", wrong dimension: (" + to_string(length) + ") expected gut got " +
                 AbstractInterface::getShapeStr(shape));
  }

  *dest = P;

}

void AbstractInterface::getOptional1DArray(float **dest,const std::string& name,int length) {

  try {
    get1DArray(dest,name,length);
  } catch (string&) {
  }

}

void AbstractInterface::get1DArray(float **dest,const std::string& name,int length) {

  *dest = nullptr;

  static vector<int64_t> shape;
  float *P = getNativeFloatArray(name,shape);
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

double AbstractInterface::getOptionalDouble(const std::string& name, double defaultValue) {

  try {
    return getDouble(name);
  } catch (string&) {
    return defaultValue;
  }

}
