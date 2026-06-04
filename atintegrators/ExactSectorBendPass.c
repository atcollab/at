#define MAGNET_PASS ExactSectorBendPass
#define INTEGRATOR_6

#include "drift_exactbend.h"
#include "kick_k1h_kn.h"
#include "curved_dipole.h"

#define INTEGRATOR(r6, num_int_steps, slength, irho, A0, B0, A, B, max_order, rad_const, diff_const, bdiff) \
    if (num_int_steps == 0) { \
        DRIFT(r6, slength, irho, bdiff); \
        FIX_LENGTH(slength); \
    } \
    else { \
        integrator(r6, num_int_steps, slength, irho, A0, B0, A, B, max_order, rad_const, diff_const, bdiff); \
    }

#define CHECK_NSTEPS \
    int ForceSplit = atGetOptionalLong(ElemData,"ForceSplit",0); check_error(); \
    int nsteps = ForceSplit ? NumIntSteps : 0; \
    for (int i=MaxOrder; i>=0; i--) { \
        if ((PolynomA[i] != 0.0) || (PolynomB[i] != 0.0)) { \
            if (NumIntSteps == 0) { \
                atError("NumIntSteps == 0 not allowed with multipoles"); check_error(); \
            } \
            nsteps = NumIntSteps; \
            break; \
        } \
    } \
    NumIntSteps=nsteps;

#include "magnet_template.h"
