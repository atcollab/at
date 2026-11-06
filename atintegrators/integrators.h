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

#if defined(INTEGRATOR_PREFIX) || defined(INTEGRATOR_SUFFIX)
#define PARTICLE_INTEGRATOR
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

#ifndef PARTICLE_BLOCK
#define PARTICLE_BLOCK 64
#endif

/* Inside the integrator, the coordinates of a block of particles are stored
   transposed (structure of arrays): coordinate k of particle c is
   rt[k*PARTICLE_BLOCK + c]. Each particle is copied into a local r6[6] array
   so that the DRIFT and KICK kernels are unchanged. Once inlined, r6 is kept
   in registers and the loop over particles accesses contiguous memory, which
   allows its vectorisation. */

/* Explicitly unrolled, so that r6 is kept in registers even at -O2 */
#define LOAD_PARTICLE(rt, c) \
    double r6[6]; \
    r6[0] = rt[c]; \
    r6[1] = rt[PARTICLE_BLOCK + c]; \
    r6[2] = rt[2*PARTICLE_BLOCK + c]; \
    r6[3] = rt[3*PARTICLE_BLOCK + c]; \
    r6[4] = rt[4*PARTICLE_BLOCK + c]; \
    r6[5] = rt[5*PARTICLE_BLOCK + c];

#define STORE_PARTICLE(rt, c) \
    rt[c] = r6[0]; \
    rt[PARTICLE_BLOCK + c] = r6[1]; \
    rt[2*PARTICLE_BLOCK + c] = r6[2]; \
    rt[3*PARTICLE_BLOCK + c] = r6[3]; \
    rt[4*PARTICLE_BLOCK + c] = r6[4]; \
    rt[5*PARTICLE_BLOCK + c] = r6[5];

/* No test of lost particles: they are processed with NaN coordinates */
#define PARTICLE_LOOP(rt, nb, ...) \
    for (int c = 0; c < nb; c++) { \
        LOAD_PARTICLE(rt, c) \
        __VA_ARGS__; \
        STORE_PARTICLE(rt, c) \
    }

#define SINGLE_PARTICLE(rt, nb, ...) __VA_ARGS__

#if defined(INTEGRATOR_4)

#define DEFAULT_NUM_INT_STEPS 10

#define INTEGRATOR_STEPS(sl) \
    double ID1 = DRIFT1 * sl; \
    double ID2 = DRIFT2 * sl; \
    double IK1 = KICK1 * sl; \
    double IK2 = KICK2 * sl;

#define integrator_slice(STEP, rt, nb, irho, A, B, max_order, rad_const, diff_const, bdiff) \
        STEP(rt, nb, DRIFT(r6, ID1, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK1, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID2, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK2, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID2, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK1, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID1, irho, bdiff));

#elif defined(INTEGRATOR_6)

#define DEFAULT_NUM_INT_STEPS 10

#define INTEGRATOR_STEPS(sl) \
    double ID1 = YD1 * sl; \
    double ID2 = YD2 * sl; \
    double ID3 = YD3 * sl; \
    double ID4 = YD4 * sl; \
    double IK1 = YK1 * sl; \
    double IK2 = YK2 * sl; \
    double IK3 = YK3 * sl; \
    double IK4 = YK4 * sl;

#define integrator_slice(STEP, rt, nb, irho, A, B, max_order, rad_const, diff_const, bdiff) \
        STEP(rt, nb, DRIFT(r6, ID1, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK1, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID2, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK2, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID3, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK3, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID4, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK4, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID4, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK3, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID3, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK2, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID2, irho, bdiff)); \
        STEP(rt, nb, KICK(r6, A, B, max_order, IK1, irho, rad_const, diff_const, bdiff)); \
        STEP(rt, nb, DRIFT(r6, ID1, irho, bdiff));

#endif /*INTEGRATOR_4*/

/* Particle by particle: lost particles are skipped, so that they do not
   contaminate the diffusion matrix nor consume random numbers */
#define integrator_particle(rt, nb, alive, num_int_steps, slength, irho, A, B, max_order, rad_const, diff_const, bdiff) \
    for (int c = 0; c < nb; c++) { \
        if (alive[c]) { \
            LOAD_PARTICLE(rt, c) \
            for (int m = 0; m < num_int_steps; m++) { /* Loop over slices */ \
                INTEGRATOR_PREFIX; \
                integrator_slice(SINGLE_PARTICLE, rt, nb, irho, A, B, max_order, rad_const, diff_const, bdiff) \
                INTEGRATOR_SUFFIX; \
            } \
            STORE_PARTICLE(rt, c) \
        } \
    }

/* Slice by slice on a block of particles: no diffusion matrix */
#define integrator_block_n(rt, nb, num_int_steps, slength, irho, A, B, max_order, rad_const, diff_const) \
    for (int m = 0; m < num_int_steps; m++) { /* Loop over slices */ \
        integrator_slice(PARTICLE_LOOP, rt, nb, irho, A, B, max_order, rad_const, diff_const, NULL) \
    }

/* Dispatch on max_order: with a constant max_order, the multipole loop in the
   inlined KICK is fully unrolled, which allows the vectorisation of the loop
   over particles. Higher orders use the generic, non-vectorised, loop. */
#define integrator_block(rt, nb, num_int_steps, slength, irho, A, B, max_order, rad_const, diff_const) \
    switch (max_order) { \
    case 0: \
        integrator_block_n(rt, nb, num_int_steps, slength, irho, A, B, 0, rad_const, diff_const) \
        break; \
    case 1: \
        integrator_block_n(rt, nb, num_int_steps, slength, irho, A, B, 1, rad_const, diff_const) \
        break; \
    case 2: \
        integrator_block_n(rt, nb, num_int_steps, slength, irho, A, B, 2, rad_const, diff_const) \
        break; \
    case 3: \
        integrator_block_n(rt, nb, num_int_steps, slength, irho, A, B, 3, rad_const, diff_const) \
        break; \
    default: \
        integrator_block_n(rt, nb, num_int_steps, slength, irho, A, B, max_order, rad_const, diff_const) \
    }

#ifdef PARTICLE_INTEGRATOR

#define integrator(rt, nb, alive, num_int_steps, slength, irho, A, B, max_order, rad_const, diff_const, bdiff) \
    integrator_particle(rt, nb, alive, num_int_steps, slength, irho, A, B, max_order, rad_const, diff_const, bdiff)

#else

#define integrator(rt, nb, alive, num_int_steps, slength, irho, A, B, max_order, rad_const, diff_const, bdiff) \
    if ((nb == 1) || bdiff) { \
        integrator_particle(rt, nb, alive, num_int_steps, slength, irho, A, B, max_order, rad_const, diff_const, bdiff) \
    } \
    else { \
        integrator_block(rt, nb, num_int_steps, slength, irho, A, B, max_order, rad_const, diff_const) \
    }

#endif /*PARTICLE_INTEGRATOR*/
