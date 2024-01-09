#include "SymplecticIntegrator.h"
#include "../atintegrators/atconstants.h"
#include <inttypes.h>
using namespace std;

SymplecticIntegrator::SymplecticIntegrator(int type) {

  switch (type) {

    case 1: // Euler
      allocate(1);
      c[0] = 1.0;
      d[0] = 1.0;
      break;

    case 2: // Verlet
      allocate(2);
      c[0] = 0.0; c[1]=1.0;
      d[0] = 0.5; d[1]=0.5;
      break;

    case 3: // Ruth
      allocate(3);
      c[0] = 1.0;       c[1]=-2.0/3.0; c[2]=2.0/3.0;
      d[0] = -1.0/24.0; d[1]=3.0/4.0;  d[2]=7.0/24.0;
      break;

    case 4: // Forest/Ruth
      allocate(4);
      c[0] = DRIFT1; c[1]=DRIFT2; c[2]=DRIFT2; c[3]=DRIFT1;
      d[0] = KICK1;  d[1]=KICK2;  d[2]=KICK1;  d[3]=0.0;
      break;

    default:
      throw string("SymplecticIntegrator, type " + to_string(type) + " not supported");

  }

}

void SymplecticIntegrator::resetMethods() {
  driftMethods.clear();
  kickMethods.clear();
}

void SymplecticIntegrator::addDriftMethod(const std::string& driftMethod) {
  driftMethods.push_back(driftMethod);
}

void SymplecticIntegrator::addKickMethod(const std::string& kickMethod) {
  kickMethods.push_back(kickMethod);
}

double SymplecticIntegrator::getLastKickWeight() {
  return d[nbCoefficients-1];
}

void SymplecticIntegrator::generateCode(std::string& code) {

  if( driftMethods.empty() || (driftMethods.size() != kickMethods.size()) )
    throw string("SymplecticIntegrator::generateCode(), wrong pass methods");

  code.append("  AT_FLOAT SL = elem->SL;\n");
  code.append("  int nb = elem->NumIntSteps;\n");

  if( driftMethods.size()==1 ) {

    generateLoopCode(code, 0);

  } else {

    code.append("  switch(elem->SubType) {\n");
    for (size_t subType = 1; subType < driftMethods.size(); subType++) {
      string sTStr = to_string(subType);
      code.append("  case " + sTStr + ":{\n");
      generateLoopCode(code, subType);
      code.append("    }break;\n");
    }
    code.append("  default:{\n");
    generateLoopCode(code, 0);
    code.append("    }break;\n");
    code.append("  }\n");

  }


}

void SymplecticIntegrator::generateLoopCode(std::string& code,size_t subType) {

  code.append("    for(int m = 0; m < nb; m++) {\n");

  for(int i=0;i<nbCoefficients;i++) {

    string drMthod = driftMethods[subType];
    if( !replace(drMthod,"%STEP%","SL*"+formatFloat(&c[i]) ) )
      throw string("SymplecticIntegrator::generateLoopCode(), wrong drift method");
    string kickMthod = kickMethods[subType];
    if( !replace(kickMthod,"%STEP%","SL*"+formatFloat(&d[i]) ) )
      throw string("SymplecticIntegrator::generateLoopCode(), wrong kick method");

    // Drift
    if( c[i]!=0.0 ) {
      code.append("      ");
      code.append(drMthod);
      code.append(";\n");
    }

    // Kick
    if( d[i]!=0.0 ) {
      code.append("      ");
      code.append(kickMthod);
      code.append(";\n");
    }

  }

  code.append("    }\n");

}

SymplecticIntegrator::~SymplecticIntegrator() {
  delete[] c;
  delete[] d;
}

void SymplecticIntegrator::allocate(int nb) {
  nbCoefficients = nb;
  c = new AT_FLOAT[nbCoefficients];
  d = new AT_FLOAT[nbCoefficients];
}

bool SymplecticIntegrator::replace(std::string& str, const std::string& from, const std::string& to) {

  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;

}

std::string SymplecticIntegrator::formatFloat(AT_FLOAT *f) {
  char bStr[128];
#ifdef WIN64
  sprintf(bStr, "__longlong_as_double(0x%016I64XULL)", *((uint64_t *)f));
#else
  sprintf(bStr, "__longlong_as_double(0x%" PRIx64  "ULL)", *((uint64_t *)f));
#endif
  return string(bStr);
}
