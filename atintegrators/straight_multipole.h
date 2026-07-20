#ifndef STRAIGHT_MULTIPOLE
#define STRAIGHT_MULTIPOLE
#include "multipolefringe.h"

#ifndef MAGNET_ENTRY
#define MAGNET_ENTRY \
    /* Entry face */ \
    multipole_fringe(r6, FringeQuadEntrance, B1, A, B, max_order, fringeIntM0, fringeIntP0, 1.0);
#endif /*MAGNET_ENTRY*/

#ifndef MAGNET_EXIT
#define MAGNET_EXIT \
    /* Exit face */ \
    multipole_fringe(r6, FringeQuadExit, B1, A, B, max_order, fringeIntM0, fringeIntP0, -1.0);
#endif /*MAGNET_EXIT*/

#define MAGNET_ARGUMENTS

#define MAGNET_ITEMS \
    Elem->BendingAngle=0.0; \
    Elem->EntranceAngle=0.0; \
    Elem->ExitAngle=0.0; \
    Elem->FringeBendEntrance=0; \
    Elem->FringeBendExit=0; \
    Elem->gK_entrance=0.0; \
    Elem->gK_exit=0.0;

#if defined(MATLAB_MEX_FILE)
#define MAGNET_MEX_ITEMS \
    double BendingAngle=0.0; \
    double EntranceAngle=0.0; \
    double ExitAngle=0.0; \
    int FringeBendEntrance=0; \
    int FringeBendExit=0; \
    double gK_entrance=0.0; \
    double gK_exit=0.0;

const char *required[] = {};
const char *optional[] = {};
#define N_REQUIRED 0
#define N_OPTIONAL 0
#endif /*MATLAB_MEX_FILE*/
#endif /*STRAIGHT_MULTIPOLE*/



