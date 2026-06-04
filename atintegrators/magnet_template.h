#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "integrators.h"

#ifndef CHECK_NSTEPS
#define CHECK_NSTEPS \
    if (NumIntSteps <= 0) { \
        atError("NumIntSteps must be positive"); check_error(); \
    }
#endif

struct elem
{
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    double BendingAngle;
    double EntranceAngle;
    double ExitAngle;
    /* Optional fields */
    double Energy;
    double Scaling;
    int FringeBendEntrance;
    int FringeBendExit;
    double gK_entrance;
    double gK_exit;
    int FringeQuadEntrance;
    int FringeQuadExit;
    double *fringeIntM0;
    double *fringeIntP0;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *KickAngle;
#ifdef STRAIGHT_DIPOLE
    double X0ref;
    double RefDZ;
#endif /*STRAIGHT_DIPOLE*/
#ifdef E2_DIPOLE
    double H1;
    double H2;
#endif /*E2_DIPOLE*/
};

static void magnet(double *r, double le, double bending_angle,
    double *A, double *B,
    int max_order, int num_int_steps,
    double entrance_angle, double exit_angle,
    int FringeBendEntrance, int FringeBendExit,
    double gK_entrance, double gK_exit,
    int FringeQuadEntrance, int FringeQuadExit,
    double *fringeIntM0,  /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
    double *fringeIntP0,  /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */
    double *T1, double *T2,
    double *R1, double *R2,
    double *RApertures, double *EApertures,
    double *KickAngle, double scaling,
#if defined(STRAIGHT_DIPOLE)
    double x0ref,
    double refdz,
#endif
#if defined(E2_DIPOLE)
    double h1,
    double h2,
#endif
#if defined(RADIATION)
    double gamma0,
    double *bdiff,
#endif
#if defined(QUANTUM)
    double gamma0,
    pcg32_random_t *rng,
#endif
    int num_particles
)
{
    double irho = bending_angle / le;

    #ifdef RADIATION
    double rad_const = RAD_CONST*pow(gamma0, 3);
    double diff_const = DIF_CONST*pow(gamma0, 5);
    #else
    double *bdiff = NULL;
    #endif

    #ifdef STRAIGHT_DIPOLE
    double phi2 = 0.5 * bending_angle;
    double phi_entrance = phi2-entrance_angle;
    double phi_exit = phi2-exit_angle;
    double LR = fabs(phi2) < 1.e-10 ? le : le *sin(phi2) / phi2;
    double SL = (num_int_steps > 0) ? LR/num_int_steps : LR;
    #else
    double refdz = 0.0;
    double SL = (num_int_steps > 0) ? le/num_int_steps : le;
    #endif /*STRAIGHT_DIPOLE*/

    INTEGRATOR_STEPS(SL)
    double B1 = (max_order >= 1) ? B[1] : 0.0;
    double A0 = 0.0;
    #ifdef CURVATURE_IN_B0
    double B0 = irho;
    #else
    double B0 = 0.0;
    #endif

    if (KickAngle) { /* Convert corrector component to polynomial coefficients */
        B0 -= sin(KickAngle[0]) / le;
        A0 += sin(KickAngle[1]) / le;
    }

    #ifndef NO_OMP
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(shared)
    #endif
    for (int c = 0; c<num_particles; c++) { /* Loop over particles */
        double *r6 = r + 6*c;
        if (!atIsNaN(r6[0])) {
            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, scaling);

            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);

            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);

            /* Entry face */
            MAGNET_ENTRY

            /* Integrator */
            INTEGRATOR(r6, num_int_steps, SL, irho, A0, B0, A, B, max_order, rad_const, diff_const, bdiff)

            /* Exit face*/
            MAGNET_EXIT

            /* Check physical apertures at the exit of the magnet */
            if (RApertures) checkiflostRectangularAp(r6, RApertures);
            if (EApertures) checkiflostEllipticalAp(r6, EApertures);

            /* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);

            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, 1.0/scaling);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double *bdiff = NULL;
    if (!Elem) {
        double Length=atGetDouble(ElemData,"Length"); check_error();
        double *PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        double *PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        int MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        int NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        MAGNET_ARGUMENTS
        /*optional fields*/
        #if defined(RADIATION) || defined(QUANTUM)
        double Energy = atEnergy(Param->energy, atGetOptionalDouble(ElemData,"Energy", Param->energy)); check_error();
        #else
        double Energy=0.0;
        #endif
        double Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        int FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error(); \
        int FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        double *fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        double *fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        double *R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        double *R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        double *T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        double *T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        double *EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        double *RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        double *KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();
        CHECK_NSTEPS

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        Elem->NumIntSteps=NumIntSteps;
        MAGNET_ITEMS
        /*optional fields*/
        Elem->Energy=Energy;
        Elem->Scaling=Scaling;
        Elem->FringeQuadEntrance=FringeQuadEntrance; \
        Elem->FringeQuadExit=FringeQuadExit; \
        Elem->fringeIntM0=fringeIntM0; \
        Elem->fringeIntP0=fringeIntP0; \
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->KickAngle=KickAngle;
    }
    #if defined(RADIATION) || defined(QUANTUM)
    double gamma0 = atGamma(Param->energy, Elem->Energy, Param->rest_energy); check_error();
    #ifdef DIFFUSION
    bdiff = Param->bdiff;
    #endif
    #endif
    magnet(r_in, Elem->Length, Elem->BendingAngle,
            Elem->PolynomA, Elem->PolynomB,
            Elem->MaxOrder, Elem->NumIntSteps, Elem->EntranceAngle, Elem->ExitAngle,
            Elem->FringeBendEntrance,Elem->FringeBendExit,
            Elem->gK_entrance, Elem->gK_exit,
            Elem->FringeQuadEntrance, Elem->FringeQuadExit,
            Elem->fringeIntM0, Elem->fringeIntP0,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            Elem->RApertures, Elem->EApertures,
            Elem->KickAngle, Elem->Scaling,
    #if defined(STRAIGHT_DIPOLE)
            Elem->X0ref,
            Elem->RefDZ,
    #endif
    #if defined(E2_DIPOLE)
            Elem->H1,
            Elem->H2,
    #endif
    #if defined(RADIATION)
            gamma0,
            bdiff,
    #endif
    #if defined(QUANTUM)
            gamma0,
            Param->thread_rng,
    #endif
            num_particles);
    return Elem;
}

MODULE_DEF(MAGNET_PASS)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgTxt("Second argument must be a 6 x N matrix");

        double Length=atGetDouble(ElemData,"Length"); check_error();
        double *PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        double *PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        int MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        int NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        MAGNET_ARGUMENTS
        /*optional fields*/
        double Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        int FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error(); \
        int FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        double *fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        double *fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        double *R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        double *R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        double *T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        double *T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        double *EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        double *RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        double *KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();
        CHECK_NSTEPS
        MAGNET_MEX_ITEMS

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);

        #if defined(RADIATION) || defined(QUANTUM)
        double Energy=atGetOptionalDouble(ElemData,"Energy",0.0); check_error();
        double rest_energy = 0.0;
        double charge = -1.0;
        if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);
        double gamma0 = atGamma(Energy, Energy, rest_energy);
        #endif

        magnet(r_in, Length, BendingAngle, PolynomA, PolynomB,
            MaxOrder, NumIntSteps, EntranceAngle, ExitAngle,
            FringeBendEntrance, FringeBendExit,
            gK_entrance, gK_exit,
            FringeQuadEntrance, FringeQuadExit,
            fringeIntM0, fringeIntP0,
            T1, T2, R1, R2, RApertures, EApertures,
            KickAngle, Scaling,
        #if defined(STRAIGHT_DIPOLE)
            X0ref,
            RefDZ,
        #endif
        #if defined(E2_DIPOLE)
            H1,
            H2,
        #endif
        #if defined(RADIATION)
            gamma0,
            NULL,
        #endif
        #if defined(QUANTUM)
            gamma0,
            &pcg32_global,
        #endif
            num_particles);
    } else if (nrhs == 0) {
        /* list of required fields */
        int i0 = 0;
        plhs[0] = mxCreateCellMatrix(5+N_REQUIRED, 1);
        mxSetCell(plhs[0], i0++, mxCreateString("Length"));
        mxSetCell(plhs[0], i0++, mxCreateString("PolynomA"));
        mxSetCell(plhs[0], i0++, mxCreateString("PolynomB"));
        mxSetCell(plhs[0], i0++, mxCreateString("MaxOrder"));
        mxSetCell(plhs[0], i0++, mxCreateString("NumIntSteps"));
        for (int i=0; i<N_REQUIRED; i++)
            mxSetCell(plhs[0], i0++, mxCreateString(required[i]));
        if (nlhs>1) {    /* list of optional fields */
            int i1 = 0;
            plhs[1] = mxCreateCellMatrix(12+N_OPTIONAL, 1);
            for (int i=0; i<N_OPTIONAL; i++)
                mxSetCell(plhs[1], i1++, mxCreateString(optional[i]));
            mxSetCell(plhs[1], i1++, mxCreateString("FringeQuadEntrance"));
            mxSetCell(plhs[1], i1++, mxCreateString("FringeQuadExit"));
            mxSetCell(plhs[1], i1++, mxCreateString("fringeIntM0"));
            mxSetCell(plhs[1], i1++, mxCreateString("fringeIntP0"));
            mxSetCell(plhs[1], i1++, mxCreateString("T1"));
            mxSetCell(plhs[1], i1++, mxCreateString("T2"));
            mxSetCell(plhs[1], i1++, mxCreateString("R1"));
            mxSetCell(plhs[1], i1++, mxCreateString("R2"));
            mxSetCell(plhs[1], i1++, mxCreateString("RApertures"));
            mxSetCell(plhs[1], i1++, mxCreateString("EApertures"));
            mxSetCell(plhs[1], i1++, mxCreateString("KickAngle"));
            mxSetCell(plhs[1], i1++, mxCreateString("FieldScaling"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
