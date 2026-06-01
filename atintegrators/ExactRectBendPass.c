#define MAGNET_PASS ExactRectBendPass
#define INTEGRATOR_6

#include "drift_exactstrbend.h"
#include "kick_exactkn.h"
#include "straight_dipole.h"

#define INTEGRATOR(r6, num_int_steps, slength, klength, irho, A0, B0, A, B, max_order, rad_const, diff_const, bdiff) \
    if (num_int_steps == 0) \
        drift(r6, slength, irho, bdiff); \
    else \
        integrator(r6, num_int_steps, slength, klength, irho, A0, B0, A, B, max_order, rad_const, diff_const, bdiff);

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
