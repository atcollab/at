#ifndef ATTYPES_H
#define ATTYPES_H

#ifndef OMP_PARTICLE_THRESHOLD
#define OMP_PARTICLE_THRESHOLD (10)
#endif

struct elem;

struct parameters
{
  int nturn;
  double RingLength;
  double T0;
  double energy;
  double rest_energy;
  int charge;
};

#endif /*ATTYPES_H*/