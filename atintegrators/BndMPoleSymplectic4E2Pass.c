#define MAGNET_PASS BndMPoleSymplectic4E2Pass
#define INTEGRATOR_4

#include "drift_E2.h"
#include "kick_E2.h"
#include "E2_dipole.h"

/*
 This code was modified from the original BndMPoleSymplectic4Pass.c of AT to correctly integrate the Hamiltonian in 
 the curvilinear coordinate system of the dipole and to include the second order Transport map of the fringe field. 
 New version created by Xiaobiao Huang in March 2009, in final verified version in August 2009.
 */

#include "magnet_template.h"
