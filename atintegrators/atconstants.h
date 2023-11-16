#ifndef AT_CONSTANTS_H
#define AT_CONSTANTS_H

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double TWOPI = 2.0*M_PI;

// Symplectic integrator constants
// Fourth-Order Symplectic Integration, E. Forest, R.D. Ruth

#define THIRD_ROOT_2 1.25992104989487316477
const double DRIFT1 = 1.0 / (2.0 * (2.0 - THIRD_ROOT_2));
const double DRIFT2 = 0.5 - (1.0 / (2.0 * (2.0 - THIRD_ROOT_2)));
const double KICK1 = 1.0 / (2.0 - THIRD_ROOT_2);
const double KICK2 = 1.0 - 2.0 / (2.0 - THIRD_ROOT_2);

// Speed of light
const double C0 = 2.99792458e8;

// Radiation damping, Physics of Electron Storage Ring, M. Sands (4.2)
#define __RE 2.8179403262e-15 // Classical electron radius [m]
#define __E0 510.99895e-6     // Electron rest energy [GeV]
const double CGAMMA = 4.0*M_PI*__RE / (3.0*__E0*__E0*__E0);

#endif //AT_CONSTANTS_H
