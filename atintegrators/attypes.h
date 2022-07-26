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
  double current;
  double nbunch;
  double bunch_spacing;
  double *fill_pattern;
  double *bunch_weights;
};

#endif /*ATTYPES_H*/
