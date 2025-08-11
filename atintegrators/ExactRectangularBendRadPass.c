#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "exactdrift.c"
#include "exactkickrad.c"
#include "exactbendfringe.c"
#include "exactmultipolefringe.c"

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
    double Energy;
    /* Optional fields */
    double Scaling;
    int FringeBendEntrance;
    int FringeBendExit;
    int FringeQuadEntrance;
    int FringeQuadExit;
    double gK;
    double x0ref;
    double refdz;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *KickAngle;
};

static void ExactRectangularBendRad(double *r, double le, double bending_angle,
        double *A, double *B,
        int max_order, int num_int_steps,
        double entrance_angle, double exit_angle,
        int FringeBendEntrance, int FringeBendExit,
        int FringeQuadEntrance, int FringeQuadExit,
        double gK, double x0ref, double refdz,
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double *KickAngle, double scaling, double gamma, int num_particles)
{
    double irho = bending_angle / le;
    double phi2 = 0.5 * bending_angle;
    double LR = phi2 < 1.e-10 ? le : le *sin(phi2) / phi2;
    double SL = LR/num_int_steps;
    double L1 = SL*DRIFT1;
    double L2 = SL*DRIFT2;
    double K1 = SL*KICK1;
    double K2 = SL*KICK2;
    double B0 = B[0];
    double A0 = A[0];
    double rad_const = RAD_CONST*pow(gamma, 3);
    double diff_const = DIF_CONST*pow(gamma, 5);

    if (KickAngle) {   /* Convert corrector component to polynomial coefficients */
        B[0] -= sin(KickAngle[0])/le;
        A[0] += sin(KickAngle[1])/le;
    }
    B[0] += irho;

    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
    shared(r,num_particles,R1,T1,R2,T2,RApertures,EApertures,\
    irho,gK,A,B,L1,L2,K1,K2,max_order,num_int_steps,scaling,\
    entrance_angle,exit_angle,x0ref,refdz,\
    FringeBendEntrance,FringeBendExit,FringeQuadEntrance,FringeQuadExit,\
    LR,le,phi2,rad_const,diff_const)
    for (int c = 0; c<num_particles; c++) { /* Loop over particles */
        double *r6 = r + 6*c;
        if (!atIsNaN(r6[0])) {
            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, scaling);

            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);

            /* Change to the magnet referential */
            Yrot(r6, entrance_angle);

            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);

            /* edge focus */
            if (FringeBendEntrance)
                bend_fringe(r6, irho, gK);
            if (FringeQuadEntrance)
                multipole_fringe(r6, le, A, B, max_order, 1.0, 1);
            bend_edge(r6, irho, phi2-entrance_angle);

            /* integrator */
            r6[0] += x0ref;
            for (int m = 0; m < num_int_steps; m++) { /* Loop over slices */
                exact_drift(r6, L1);
                ex_strthinkickrad(r6, A, B, max_order, 0.0, K1, rad_const, diff_const, NULL);
                exact_drift(r6, L2);
                ex_strthinkickrad(r6, A, B, max_order, 0.0, K2, rad_const, diff_const, NULL);
                exact_drift(r6, L2);
                ex_strthinkickrad(r6, A, B, max_order, 0.0, K1, rad_const, diff_const, NULL);
                exact_drift(r6, L1);
            }
            r6[0] -= x0ref;

            /* Convert absolute path length to path lengthening */
            r6[5] -= (le+refdz);

            /* edge focus */
            bend_edge(r6, irho, phi2-exit_angle);
            if (FringeQuadExit)
                multipole_fringe(r6, le, A, B, max_order, -1.0, 1);
            if (FringeBendExit)
                bend_fringe(r6, -irho, gK);

            /* Check physical apertures at the exit of the magnet */
            if (RApertures) checkiflostRectangularAp(r6, RApertures);
            if (EApertures) checkiflostEllipticalAp(r6, EApertures);

            /* Change back to the lattice referential */
            Yrot(r6, exit_angle);

            /* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);

            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, 1.0/scaling);
        }
    }
    /* Remove corrector component in polynomial coefficients */
    B[0] = B0;
    A[0] = A0;
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double gamma;
    if (!Elem) {
        double Length=atGetDouble(ElemData,"Length"); check_error();
        double *PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        double *PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        int MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        int NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        double BendingAngle=atGetOptionalDouble(ElemData,"BendingAngle", 0.0); check_error();
        double EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        double ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        double Energy=atGetOptionalDouble(ElemData,"Energy", Param->energy); check_error();
        double Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        int FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        int FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        int FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error();
        int FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        double gK=atGetOptionalDouble(ElemData,"gK", 0.0); check_error();
        double x0ref=atGetOptionalDouble(ElemData,"X0ref", 0.0); check_error();
        double refdz=atGetOptionalDouble(ElemData,"RefDZ", 0.0); check_error();
        double *R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        double *R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        double *T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        double *T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        double *EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        double *RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        double *KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();

        if (NumIntSteps <= 0) {
            atError("NumIntSteps must be positive"); check_error();
        }

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        Elem->NumIntSteps=NumIntSteps;
        Elem->BendingAngle=BendingAngle;
        Elem->EntranceAngle=EntranceAngle;
        Elem->ExitAngle=ExitAngle;
        /*optional fields*/
        Elem->Energy=Energy;
        Elem->Scaling=Scaling;
        Elem->FringeBendEntrance=FringeBendEntrance;
        Elem->FringeBendExit=FringeBendExit;
        Elem->FringeQuadEntrance=FringeQuadEntrance;
        Elem->FringeQuadExit=FringeQuadExit;
        Elem->gK=gK;
        Elem->x0ref=x0ref;
        Elem->refdz=refdz;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->KickAngle=KickAngle;
    }
    gamma = atGamma(Param->energy, Elem->Energy, Param->rest_energy);

    ExactRectangularBendRad(r_in, Elem->Length, Elem->BendingAngle, Elem->PolynomA, Elem->PolynomB,
            Elem->MaxOrder, Elem->NumIntSteps, Elem->EntranceAngle, Elem->ExitAngle,
            Elem->FringeBendEntrance,Elem->FringeBendExit,
            Elem->FringeQuadEntrance, Elem->FringeQuadExit,
            Elem->gK,Elem->x0ref,Elem->refdz,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            Elem->RApertures, Elem->EApertures,
            Elem->KickAngle, Elem->Scaling, gamma, num_particles);
    return Elem;
}

MODULE_DEF(ExactRectangularBendRadPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double rest_energy = 0.0;
        double charge = -1.0;
        double Gamma;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgTxt("Second argument must be a 6 x N matrix");

        double Length=atGetDouble(ElemData,"Length"); check_error();
        double *PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        double *PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        int MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        int NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        double BendingAngle=atGetOptionalDouble(ElemData,"BendingAngle",0.0); check_error();
        double EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        double ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        double Energy=atGetOptionalDouble(ElemData,"Energy",0.0); check_error();
        double Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        int FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        int FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        int FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error();
        int FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        double gK=atGetOptionalDouble(ElemData,"gK", 0.0); check_error();
        double x0ref=atGetOptionalDouble(ElemData,"X0ref", 0.0); check_error();
        double refdz=atGetOptionalDouble(ElemData,"RefDZ", 0.0); check_error();
        double *R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        double *R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        double *T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        double *T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        double *EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        double *RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        double *KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();

        if (NumIntSteps <= 0) {
            atError("NumIntSteps must be positive"); check_error();
        }
        if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        Gamma = atGamma(Energy, Energy, rest_energy);
        r_in = mxGetDoubles(plhs[0]);

        ExactRectangularBendRad(r_in, Length, BendingAngle, PolynomA, PolynomB,
            MaxOrder, NumIntSteps, EntranceAngle, ExitAngle,
            FringeBendEntrance, FringeBendExit,
            FringeQuadEntrance, FringeQuadExit,
            gK, x0ref, refdz,
            T1, T2, R1, R2, RApertures, EApertures,
            KickAngle, Scaling, Gamma, num_particles);
    } else if (nrhs == 0) {
        /* list of required fields */
        int i0 = 0;
        plhs[0] = mxCreateCellMatrix(8, 1);
        mxSetCell(plhs[0], i0++, mxCreateString("Length"));
        mxSetCell(plhs[0], i0++, mxCreateString("PolynomA"));
        mxSetCell(plhs[0], i0++, mxCreateString("PolynomB"));
        mxSetCell(plhs[0], i0++, mxCreateString("MaxOrder"));
        mxSetCell(plhs[0], i0++, mxCreateString("NumIntSteps"));
        mxSetCell(plhs[0], i0++, mxCreateString("BendingAngle"));
        mxSetCell(plhs[0], i0++, mxCreateString("EntranceAngle"));
        mxSetCell(plhs[0], i0++, mxCreateString("ExitAngle"));
        if (nlhs>1) {    /* list of optional fields */
            int i1 = 0;
            plhs[1] = mxCreateCellMatrix(16, 1);
            mxSetCell(plhs[1], i1++, mxCreateString("Energy"));
            mxSetCell(plhs[1], i1++, mxCreateString("FringeBendEntrance"));
            mxSetCell(plhs[1], i1++, mxCreateString("FringeBendExit"));
            mxSetCell(plhs[1], i1++, mxCreateString("FringeQuadEntrance"));
            mxSetCell(plhs[1], i1++, mxCreateString("FringeQuadExit"));
            mxSetCell(plhs[1], i1++, mxCreateString("gK"));
            mxSetCell(plhs[1], i1++, mxCreateString("X0ref"));
            mxSetCell(plhs[1], i1++, mxCreateString("RefDZ"));
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
