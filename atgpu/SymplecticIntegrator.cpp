#include "SymplecticIntegrator.h"
#include "AbstractGPU.h"
#include "../atintegrators/atconstants.h"
using namespace std;

SymplecticIntegrator::SymplecticIntegrator(int type) {
  c = nullptr;
  d = nullptr;
  nbCoefficients = 0;
  setType(type);
}

SymplecticIntegrator::~SymplecticIntegrator() {
  delete[] c;
  delete[] d;
}

void SymplecticIntegrator::setType(int type) {

  delete[] c;
  delete[] d;
  c = nullptr;
  d = nullptr;
  nbCoefficients = 0;
  this->type = type;

  // Note:
  // The integrators below have been tested on hamiltonian using quadratic approximation
  // for the transverse momenta. McLachlan/Atela integrators give good results.
  // TODO: Test on exact pass

  switch (type) {

    case 1: // Euler
      allocate(1);
      c[0] = 1.0;
      d[0] = 1.0;
      break;

    case 2: // Optimal 2nd order from "The Accuracy of Symplectic Integrators", R. Mclachlan P Atela, 1991
      allocate(2);
      c[0] = 1.0 - 1.0 / sqrt(2.0); c[1]=1.0 / sqrt(2.0);
      d[0] = 1.0 - 1.0 / sqrt(2.0); d[1]=1.0 / sqrt(2.0);
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

    case 5: // Optimal 4th order from "The Accuracy of Symplectic Integrators", R. Mclachlan P Atela, 1991
      allocate(4);
      c[0]= 0.1288461583653843; c[1]= 0.4415830236164665; c[2]=-0.0857820194129737; c[3]= 0.5153528374311229;
      d[0]= 0.3340036032863214; d[1]= 0.7563200005156683; d[2]=-0.2248198030794208; d[3]= 0.1344961992774311;
      break;

    case 6: // H. Yoshida, "Construction of higher order symplectic integrators", PHYSICS LETTERS A 12/11/1990
      // 6th order symplectic integrator constants (solution A)
      //w1=-0.117767998417887e1
      //w2= 0.235573213359357e0
      //w3= 0.784513610477560e0
      //w0= 1.0 - 2.0*(w1+w2+w3)
      //# Drift coef
      //print(f"c[0]= {0.5*w3:.16f}; c[1]= {0.5*(w3+w2):.16f}; c[2]={0.5*(w2+w1):.16f}; c[3]= {0.5*(w1+w0):.16f};")
      //print(f"c[4]= {0.5*(w1+w0):.16f}; c[5]={0.5*(w2+w1):.16f}; c[6]= {0.5*(w3+w2):.16f}; c[7]= {0.5*w3:.16f};")
      //# Kick coef
      //print(f"d[0]= {w3:.16f}; d[1]= {w2:.16f}; d[2]={w1:.16f}; d[3]= {w0:.16f};")
      //print(f"d[4]={w1:.16f}; d[5]= {w2:.16f}; d[6]= {w3:.16f}; d[7]= 0.0;")
      allocate(8);
      c[0]= 0.3922568052387800; c[1]= 0.5100434119184585; c[2]=-0.4710533854097566; c[3]= 0.0687531682525181;
      c[4]= 0.0687531682525181; c[5]=-0.4710533854097566; c[6]= 0.5100434119184585; c[7]= 0.3922568052387800;
      d[0]= 0.7845136104775600; d[1]= 0.2355732133593570; d[2]=-1.1776799841788701; d[3]= 1.3151863206839063;
      d[4]=-1.1776799841788701; d[5]= 0.2355732133593570; d[6]= 0.7845136104775600; d[7]= 0.0;
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

  if( driftMethods.empty() )
    return;
  if(driftMethods.size() != kickMethods.size())
    throw string("SymplecticIntegrator::generateCode(), wrong pass methods");

  code.append("  AT_FLOAT SL = elem->SL;\n");
  code.append("  int nb = elem->mpole.NumIntSteps;\n");

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
    if( !replace(drMthod,"%STEP%","SL*"+AbstractGPU::getInstance()->formatFloat(&c[i]) ) )
      throw string("SymplecticIntegrator::generateLoopCode(), wrong drift method");
    string kickMthod = kickMethods[subType];
    if( !replace(kickMthod,"%STEP%","SL*"+AbstractGPU::getInstance()->formatFloat(&d[i]) ) )
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

void SymplecticIntegrator::allocate(int nb) {
  nbCoefficients = nb;
  c = new double[nbCoefficients];
  d = new double[nbCoefficients];
}

bool SymplecticIntegrator::replace(std::string& str, const std::string& from, const std::string& to) {

  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;

}

