#ifndef STRAIGHT_DIPOLE
#define STRAIGHT_DIPOLE
#include "bendfringe.h"
#include "multipolefringe.h"

#ifndef MAGNET_ENTRY
#define MAGNET_ENTRY \
    /* Change to the magnet referential */ \
    Yrot(r6, entrance_angle, bdiff); \
    \
    /* Entry face */ \
    r6[0] += x0ref; \
    if (FringeBendEntrance) \
        bend_fringe(r6, irho, gK_entrance); \
    multipole_fringe(r6, FringeQuadEntrance, B1, A, B, max_order, fringeIntM0, fringeIntP0, 1.0); \
    if (phi_entrance != 0.0) { \
        if (B1 != 0.0 && FringeBendEntrance) quad_wedge(r6, -B1 * phi_entrance); \
        bend_wedge(r6, irho, phi_entrance, bdiff); \
    }
#endif /*MAGNET_ENTRY*/

#ifndef MAGNET_EXIT
#define MAGNET_EXIT \
    /* Exit face */ \
    if (phi_exit != 0.0) { \
        bend_wedge(r6, irho, phi_exit, bdiff); \
        if (B1 != 0.0 && FringeQuadExit) quad_wedge(r6, -B1 * phi_exit); \
    } \
    multipole_fringe(r6, FringeQuadExit, B1, A, B, max_order, fringeIntM0, fringeIntP0, -1.0); \
    if (FringeBendExit) \
        bend_fringe(r6, -irho, gK_exit); \
    r6[0] -= x0ref; \
    \
    /* Change back to the lattice referential */ \
    Yrot(r6, exit_angle, bdiff);
#endif /*MAGNET_EXIT*/

#define MAGNET_ARGUMENTS \
    double BendingAngle=atGetOptionalDouble(ElemData,"BendingAngle", 0.0); check_error(); \
    double EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error(); \
    double ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error(); \
    int FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",4); check_error(); \
    int FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",4); check_error(); \
    double FullGap=atGetOptionalDouble(ElemData,"FullGap",0.0); check_error(); \
    double FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0.0); check_error(); \
    double FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0.0); check_error(); \
    double X0ref=atGetOptionalDouble(ElemData,"X0ref", 0.0); check_error(); \
    double RefDZ=atGetOptionalDouble(ElemData,"RefDZ", 0.0); check_error();

#define MAGNET_ITEMS \
    Elem->BendingAngle=BendingAngle; \
    Elem->EntranceAngle=EntranceAngle; \
    Elem->ExitAngle=ExitAngle; \
    Elem->FringeBendEntrance=FringeBendEntrance; \
    Elem->FringeBendExit=FringeBendExit; \
    Elem->gK_entrance=FullGap*FringeInt1; \
    Elem->gK_exit=FullGap*FringeInt2; \
    Elem->X0ref=X0ref; \
    Elem->RefDZ=RefDZ;

#ifdef MATLAB_MEX_FILE
#define MAGNET_MEX_ITEMS \
    double gK_entrance=FullGap*FringeInt1; \
    double gK_exit=FullGap*FringeInt2;

const char *required[] = {"BendingAngle", "EntranceAngle", "ExitAngle"};
const char *optional[] = {"FringeBendEntrance", "FringeBendExit", "FullGap", "FringeInt1", "FringeInt2", "X0ref", "RefDZ"};
#define N_REQUIRED 3
#define N_OPTIONAL 7
#endif /*MATLAB_MEX_FILE*/
#endif /*STRAIGHT_DIPOLE*/



