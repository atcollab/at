#define MAGNET_PASS ExactSectorBendQuantPass
#define INTEGRATOR_6
#define QUANTUM
#define NO_OMP  /* because of problems with random generator and OpenMP */

#include "drift_exactbend.h"
#include "kick_k1h_kn.h"
#include "curved_dipole.h"

#include "magnet_template.h"
