#ifndef BENT_DIPOLE
#define BENT_DIPOLE
#include "bendfringe.h"
#include "multipolefringe.h"

#ifndef DEFAULT_BEND_FRINGE
#define DEFAULT_BEND_FRINGE 4
#endif

#ifndef MAGNET_ENTRY
#define MAGNET_ENTRY \
    /* Entry face */ \
    if (FringeBendEntrance == 4) { \
        Yrot(r6, entrance_angle, bdiff); \
        bend_fringe(r6, irho, gK_entrance); \
        multipole_fringe(r6, FringeQuadEntrance, B1, A, B, max_order, fringeIntM0, fringeIntP0, 1.0); \
        if (entrance_angle != 0.0) { \
            if (B1 != 0.0 && FringeQuadEntrance) quad_wedge(r6, -B1 * entrance_angle); \
            bend_wedge(r6, irho, -entrance_angle, bdiff); \
        } \
    } \
    else if (FringeBendEntrance > 0) { \
        bend_linear_fringe(r6, irho, entrance_angle, gK_entrance, FringeBendEntrance, 1.0, bdiff); \
        multipole_fringe(r6, FringeQuadEntrance, B1, A, B, max_order, fringeIntM0, fringeIntP0, 1.0); \
    }
#endif /*MAGNET_ENTRY*/

#ifndef MAGNET_EXIT
#define MAGNET_EXIT \
    /* Exit face */ \
    if (FringeBendExit == 4) { \
        if (exit_angle != 0.0) { \
            bend_wedge(r6, irho, -exit_angle, bdiff); \
            if (B1 != 0.0 && FringeQuadExit) quad_wedge(r6, -B1 * exit_angle); \
        } \
        multipole_fringe(r6, FringeQuadExit, B1, A, B, max_order, fringeIntM0, fringeIntP0, -1.0); \
        bend_fringe(r6, -irho, gK_exit); \
        Yrot(r6, exit_angle, bdiff); \
    } \
    else if (FringeBendExit > 0) { \
        multipole_fringe(r6, FringeQuadExit, B1, A, B, max_order, fringeIntM0, fringeIntP0, -1.0); \
        bend_linear_fringe(r6, irho, exit_angle, gK_exit, FringeBendExit, -1.0, bdiff); \
    }
#endif /*MAGNET_EXIT*/

#define MAGNET_ARGUMENTS \
    double BendingAngle=atGetOptionalDouble(ElemData,"BendingAngle", 0.0); check_error(); \
    double EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error(); \
    double ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error(); \
    int FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",DEFAULT_BEND_FRINGE); check_error(); \
    int FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",DEFAULT_BEND_FRINGE); check_error(); \
    double FullGap=atGetOptionalDouble(ElemData,"FullGap",0.0); check_error(); \
    double FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0.0); check_error(); \
    double FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0.0); check_error();

#define MAGNET_ITEMS \
    Elem->BendingAngle=BendingAngle; \
    Elem->EntranceAngle=EntranceAngle; \
    Elem->ExitAngle=ExitAngle; \
    Elem->FringeBendEntrance=FringeBendEntrance; \
    Elem->FringeBendExit=FringeBendExit; \
    Elem->gK_entrance=FullGap*FringeInt1; \
    Elem->gK_exit=FullGap*FringeInt2;

#ifdef MATLAB_MEX_FILE
#define MAGNET_MEX_ITEMS \
    double gK_entrance=FullGap*FringeInt1; \
    double gK_exit=FullGap*FringeInt2;

const char *required[] = {"BendingAngle", "EntranceAngle", "ExitAngle"};
const char *optional[] = {"FringeBendEntrance", "FringeBendExit", "FullGap", "FringeInt1", "FringeInt2"};
#define N_REQUIRED 3
#define N_OPTIONAL 5
#endif /*MATLAB_MEX_FILE*/
#endif /*BENT_DIPOLE*/
