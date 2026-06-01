#ifndef E2_DIPOLE
#define E2_DIPOLE
#include <math.h>
#include "bendfringe.h"
#include "multipolefringe.h"

static void edge_fringe2A(double* r, double h, double edge_angle, double gK, double h1, double K1)
{   /* Entrance Fringe field transport map to second order in dipoles with fringe field */
    double ca = cos(edge_angle);
    double sa = sin(edge_angle);
    double dpsi = h*gK*(1.0+sa*sa)/ca; /* /(1+r[4]); */
    double psi_bar = edge_angle-dpsi;
    double tpsi=sa/ca, tpsib=tan(psi_bar);
    double fx = h*tpsi;
    double fy = h*tpsib;
    double spsi=1.0/ca; /* spsib=1.0/cos(psi_bar) */
    double T111,T234,T414,   T212,T313,  T133,T423,T211,T233,T413;

    double r0=r[0],r2=r[2],r1=r[1];
    T111 = -0.5*h*tpsi*tpsi;
    /*  T234=  -0.5*h*tpsi*tpsib;    */
    T234=  -0.5*h*tpsi*tpsi;
    T414=T234;
    T212 = -T111;
    T313 = -T234;
    T133 =  0.5*h*spsi*spsi;  T423=-T133;
    T211 =  0.5*h*h1*spsi*spsi*spsi + K1*tpsi;
    T233 = -0.5*h*h1*spsi*spsi*spsi -K1*tpsi+0.5*h*h*tpsi*(tpsib*tpsib+spsi*spsi);
    T413 = -0.5*h*h1*spsi*spsi*spsi -K1*tpsi; /*-0.5*h*h*tpsi*(spsi*spsi+tpsib*tpsib);*/

    r[0] += T111*r[0]*r[0]+T133*r[2]*r[2];
    r[1] += r0*fx + 2*T212*r0*r[1]+2*T234*r[2]*r[3]+T211*r0*r0+T233*r[2]*r[2] ;
    r[2] += 2*T313*r0*r[2];
    r[3] += -r2*fy + 2*T414*r0*r[3]+2*T413*r0*r2+2*T423*r1*r2 ;

}

static void edge_fringe2B(double* r, double h, double edge_angle, double gK, double h2,double K1)
{   /* Exit Fringe field transport map to second order in dipoles with fringe field */
    double ca = cos(edge_angle);
    double sa = sin(edge_angle);
    double dpsi = h*gK*(1.0+sa*sa)/ca; /* /(1+r[4]); */
    double psi_bar = edge_angle-dpsi;
    double tpsi=sa/ca, tpsib=tan(psi_bar);
    double fx = h*tpsi;
    double fy = h*tpsib;
    double spsi=1.0/ca; /* spsib=1.0/cos(psi_bar) */
    double T111,T234,T414,   T212,T313,  T133,T423,T211,T233,T413;

    double r0=r[0],r2=r[2],r1=r[1];
    T111 = 0.5*h*tpsi*tpsi;
    /*  T234=  0.5*h*tpsi*tpsib;    */
    T234=  0.5*h*tpsi*tpsi;
    T414=T234;
    T212 = -T111;
    T313 = -T234;
    T133 = -0.5*h*spsi*spsi;  T423=-T133;
    T211 =  0.5*h*h2*spsi*spsi*spsi +K1*tpsi-0.5*h*h*tpsi*tpsi*tpsi;
    T233 = -0.5*h*h2*spsi*spsi*spsi -K1*tpsi-0.5*h*h*tpsi*tpsib*tpsib;
    T413 = -0.5*h*h2*spsi*spsi*spsi -K1*tpsi+0.5*h*h*tpsi*(spsi*spsi);

    r[0] += T111*r[0]*r[0]+T133*r[2]*r[2];
    r[1] += r0*fx + 2*T212*r0*r[1]+2*T234*r[2]*r[3]+T211*r0*r0+T233*r[2]*r[2] ;
    r[2] += 2*T313*r0*r[2];
    r[3] += -r2*fy + 2*T414*r0*r[3]+2*T413*r0*r2+2*T423*r1*r2 ;

}

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
        edge_fringe2A(r6, irho, entrance_angle, gK_entrance, h1, B1); \
        multipole_fringe(r6, FringeQuadEntrance, B1, A, B, max_order, fringeIntM0, fringeIntP0, 1.0); \
    }

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
        edge_fringe2B(r6, irho, exit_angle, gK_exit, h2, B1); \
    }

#define MAGNET_ARGUMENTS \
    double BendingAngle=atGetOptionalDouble(ElemData,"BendingAngle", 0.0); check_error(); \
    double EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error(); \
    double ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error(); \
    int FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance", 1); check_error(); \
    int FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit", 1); check_error(); \
    double FullGap=atGetOptionalDouble(ElemData,"FullGap",0.0); check_error(); \
    double FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0.0); check_error(); \
    double FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0.0); check_error(); \
    double H1=atGetOptionalDouble(ElemData,"H1",0.0); check_error(); \
    double H2=atGetOptionalDouble(ElemData,"H2",0.0); check_error();

#define MAGNET_ITEMS \
    Elem->BendingAngle=BendingAngle; \
    Elem->EntranceAngle=EntranceAngle; \
    Elem->ExitAngle=ExitAngle; \
    Elem->FringeBendEntrance=FringeBendEntrance; \
    Elem->FringeBendExit=FringeBendExit; \
    Elem->gK_entrance=FullGap*FringeInt1; \
    Elem->gK_exit=FullGap*FringeInt2; \
    Elem->H1 = H1; \
    Elem->H2 = H2;

#ifdef MATLAB_MEX_FILE
#define MAGNET_MEX_ITEMS \
    double gK_entrance=FullGap*FringeInt1; \
    double gK_exit=FullGap*FringeInt2;

const char *required[] = {"BendingAngle", "EntranceAngle", "ExitAngle"};
const char *optional[] = {"FringeBendEntrance", "FringeBendExit", "FullGap", "FringeInt1", "FringeInt2", "H1", "H2"};
#define N_REQUIRED 3
#define N_OPTIONAL 7
#endif /*MATLAB_MEX_FILE*/

#endif /*E2_DIPOLE*/
