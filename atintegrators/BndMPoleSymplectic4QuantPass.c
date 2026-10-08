#define MAGNET_PASS BndMPoleSymplectic4QuantPass
#define DEFAULT_BEND_FRINGE 1
#define INTEGRATOR_4
#define QUANTUM
#define NO_OMP  /* because of problems with random generator and OpenMP */

#include "drift_expanded.h"
#include "kick_h_k0h_k1h_kn.h"
#include "curved_dipole.h"

#include "magnet_template.h"
