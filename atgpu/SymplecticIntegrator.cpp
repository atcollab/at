#include "SymplecticIntegrator.h"
#include "../atintegrators/atconstants.h"

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

SymplecticIntegrator::~SymplecticIntegrator() {
  delete[] c;
  delete[] d;
}

void SymplecticIntegrator::allocate(int nb) {
  nbCoefficients = nb;
  c = new AT_FLOAT[nbCoefficients];
  d = new AT_FLOAT[nbCoefficients];
}
