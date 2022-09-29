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
  double s_coord;
  double energy;
  double rest_energy;
  double charge;
  double beam_current;
  int nbunch;
  double *bunch_spos;
  double *bunch_currents;
};

#endif /*ATTYPES_H*/
