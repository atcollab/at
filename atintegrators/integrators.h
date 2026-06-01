#ifdef QUANTUM
#include "quantum_diffusion.h"
#endif /*QUANTUM*/

#define YD1 3.922568052387799819591407413100e-01
#define YD2 5.100434119184584780271052295575e-01
#define YD3 -4.710533854097565531482416645304e-01
#define YD4 6.875316825251809316199569366290e-02

#define YK1 7.845136104775599639182814826199e-01
#define YK2 2.355732133593569921359289764951e-01
#define YK3 -1.177679984178870098432412305556e+00
#define YK4 1.315186320683906284756403692882e+00

#ifndef DRIFT
#define DRIFT drift
#endif

#ifndef KICK
#define KICK kick
#endif

#ifdef RADIATION
#define KICK_(r6, A0, B0, A, B, max_order, length, irho, rad_const, diff_const, bdiff) \
    KICK(r6, A0, B0, A, B, max_order, length, irho, rad_const, diff_const, bdiff)
#else
#define KICK_(r6, A0, B0, A, B, max_order, length, irho, rad_const, diff_const, bdiff) \
    KICK(r6, A0, B0, A, B, max_order, length, irho)
#endif

#ifndef INTEGRATOR_PREFIX
#define INTEGRATOR_PREFIX
#endif

#ifndef INTEGRATOR_SUFFIX
#define INTEGRATOR_SUFFIX
#endif

#ifndef INTEGRATOR
#define INTEGRATOR integrator
#endif

#ifdef ABSOLUTE_PATH_LENGTH
#define FIX_LENGTH r6[5] -= (le+refdz);
#else
#define FIX_LENGTH
#endif

#if defined(INTEGRATOR_4)
#define integrator(r6, num_int_steps, slength, klength, irho, A0, B0, A, B, max_order, rad_const, diff_const, bdiff) \
    for (int m = 0; m < num_int_steps; m++) { /* Loop over slices */ \
        INTEGRATOR_PREFIX \
        DRIFT(r6, DRIFT1 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, KICK1 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, DRIFT2 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, KICK2 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, DRIFT2 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, KICK1 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, DRIFT1 * slength, irho, bdiff); \
        INTEGRATOR_SUFFIX \
    } \
    FIX_LENGTH
#elif defined(INTEGRATOR_6)
#define integrator(r6, num_int_steps, slength, klength, irho, A0, B0, A, B, max_order, rad_const, diff_const, bdiff) \
    for (int m = 0; m < num_int_steps; m++) { /* Loop over slices */ \
        INTEGRATOR_PREFIX \
        DRIFT(r6, YD1 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, YK1 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, YD2 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, YK2 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, YD3 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, YK3 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, YD4 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, YK4 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, YD4 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, YK3 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, YD3 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, YK2 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, YD2 * slength, irho, bdiff); \
        KICK_(r6, A0, B0, A, B, max_order, YK1 * klength, irho, rad_const, diff_const, bdiff); \
        DRIFT(r6, YD1 * slength, irho, bdiff); \
        INTEGRATOR_SUFFIX \
    } \
    FIX_LENGTH
#endif /*INTEGRATOR_4*/
