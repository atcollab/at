#define MAGNET_PASS ExactMultipoleQuantPass
#define INTEGRATOR_6
#define QUANTUM
#define NO_OMP  /* because of problems with random generator and OpenMP */

#include "drift_exact.h"
#include "kick_exactkn.h"  /* kick */
#include "straight_multipole.h"

#include "magnet_template.h"
