#define MAGNET_PASS StrMPoleSymplectic6QuantPass
#define INTEGRATOR_6
#define QUANTUM
#define NO_OMP  /* because of problems with random generator and OpenMP */

#include "drift_expanded.h"
#include "kick_kn.h"
#include "straight_multipole.h"

#include "magnet_template.h"
