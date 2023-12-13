#include "SymplecticIntegrator.h"
#include "../atintegrators/atconstants.h"
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

  }

}

void SymplecticIntegrator::generateCode(string& code,const string& count,const string& driftMethod,
                                        const string& kickMethod, const string& exKickParam) {

  code.append("    for(int m = 0; m < ");
  code.append(count);
  code.append("; m++) {\n");

  for(int i=0;i<nbCoefficients;i++) {

    string paramDriftList = "(r6,elem->NormD[" + to_string(i) + "] * p_norm,p_norm);";
    string paramKickList = "(r6,elem->PolynomA,elem->PolynomB,elem->NormK[" + to_string(i) + "],elem->MaxOrder";
    paramKickList += exKickParam.empty()?");":","+exKickParam+");";

    // Drift
    if( c[i]!=0.0 ) {
      code.append("      ");
      code.append(driftMethod);
      code.append(paramDriftList);
      code.append("\n");
    }

    // Kick
    if( d[i]!=0.0 ) {
      code.append("      ");
      code.append(kickMethod);
      code.append(paramKickList);
      code.append("\n");
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
