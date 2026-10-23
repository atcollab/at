#define MAGNET_PASS BndMPoleSymplectic4E2RadPass
#define INTEGRATOR_4
#define RADIATION

#include "drift_E2.h"
#include "kick_E2.h"
#include "E2_dipole.h"

/*
 This code was modified from the original BndMPoleSymplectic4RadPass.c of AT to correctly integrate the Hamiltonian in
 the curvilinear coordinate system of the dipole and to include the second order Transport map of the fringe field. Also
 modified is the field Bx, By to include the curvature effect.
 New version created by Xiaobiao Huang on 08/13/2009.
 Last modified on 8/26/2009
 */

#include "magnet_template.h"
