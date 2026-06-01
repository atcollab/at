#define MAGNET_PASS StrMPoleSymplectic4QuantPass
#define INTEGRATOR_4
#define QUANTUM
#define NO_OMP  /* because of problems with random generator and OpenMP */

#include "drift_fast.h"    /* drift */
#include "kick_kn.h"  /* kick */
#include "straight_multipole.h"

#include "magnet_template.h"
