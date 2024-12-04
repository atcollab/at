#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "driftkick.c"  /* strthinkick.c */
#include "exactbend.c"
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
    /* Optional fields */
    double Scaling;
    int FringeBendEntrance;
    int FringeBendExit;
    int FringeQuadEntrance;
    int FringeQuadExit;
    double gK;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *KickAngle;
};

static void ExactSectorBend(double *r, double le, double bending_angle,
        double *A, double *B,
        int max_order, int num_int_steps,
        double entrance_angle, double exit_angle,
        int FringeBendEntrance, int FringeBendExit,
        int FringeQuadEntrance, int FringeQuadExit,
        double gK,
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double *KickAngle, double scaling, int num_particles)
{
    double irho = bending_angle / le;
    double SL = le/num_int_steps;
    double L1 = SL*DRIFT1;
    double L2 = SL*DRIFT2;
    double K1 = SL*KICK1;
    double K2 = SL*KICK2;
    double B0 = B[0];
    double A0 = A[0];

    if (KickAngle) {   /* Convert corrector component to polynomial coefficients */
        B[0] -= sin(KickAngle[0])/le;
        A[0] += sin(KickAngle[1])/le;
    }

    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
    shared(r,num_particles,R1,T1,R2,T2,RApertures,EApertures,\
    irho,gK,A,B,L1,L2,K1,K2,max_order,num_int_steps,scaling,\
    entrance_angle,exit_angle,\
    FringeBendEntrance,FringeBendExit,FringeQuadEntrance,FringeQuadExit,le)
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

            Yrot(r6, entrance_angle);
            if (FringeBendEntrance)
                bend_fringe(r6, irho, gK);
            if (FringeQuadEntrance)
                multipole_fringe(r6, le, A, B, max_order, 1.0, 1);
            bend_edge(r6, irho, -entrance_angle);

            if (num_int_steps == 0) {
                exact_bend(r6, irho, le);
            }
            else {
                for (int m = 0; m < num_int_steps; m++) { /* Loop over slices */
                    exact_bend(r6, irho, L1);
                    strthinkick(r6, A, B, K1, max_order);
                    exact_bend(r6, irho, L2);
                    strthinkick(r6, A, B, K2, max_order);
                    exact_bend(r6, irho, L2);
                    strthinkick(r6, A, B, K1, max_order);
                    exact_bend(r6, irho, L1);
                }
            }

            /* Convert absolute path length to path lengthening */
            r6[5] -= le;

            /* edge focus */
            bend_edge(r6, irho, -exit_angle);
            if (FringeQuadExit)
                multipole_fringe(r6, le, A, B, max_order, -1.0, 1);
            if (FringeBendExit)
                bend_fringe(r6, -irho, gK);
            Yrot(r6, exit_angle);

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
    /* Remove corrector component in polynomial coefficients */
    B[0] = B0;
    A[0] = A0;
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
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
        double Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        int FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        int FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        int FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error();
        int FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        double gK=atGetOptionalDouble(ElemData,"gK", 0.0); check_error();
        double *R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        double *R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        double *T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        double *T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        double *EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        double *RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        double *KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();

        if (NumIntSteps == 0) {
            for (int i=MaxOrder; i>=0; i--) {
                if ((PolynomA[i] != 0.0) || (PolynomB[i] != 0.0)) {
                    atError("NumIntSteps == 0 not allowed with multipoles"); check_error();
                }
            }
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
        Elem->Scaling=Scaling;
        Elem->FringeBendEntrance=FringeBendEntrance;
        Elem->FringeBendExit=FringeBendExit;
        Elem->FringeQuadEntrance=FringeQuadEntrance;
        Elem->FringeQuadExit=FringeQuadExit;
        Elem->gK=gK;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->KickAngle=KickAngle;
    }
    ExactSectorBend(r_in, Elem->Length, Elem->BendingAngle,
            Elem->PolynomA, Elem->PolynomB,
            Elem->MaxOrder, Elem->NumIntSteps, Elem->EntranceAngle, Elem->ExitAngle,
            Elem->FringeBendEntrance,Elem->FringeBendExit,
            Elem->FringeQuadEntrance, Elem->FringeQuadExit,
            Elem->gK,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            Elem->RApertures, Elem->EApertures,
            Elem->KickAngle, Elem->Scaling, num_particles);
    return Elem;
}

MODULE_DEF(ExactSectorBendPass)        /* Dummy module initialisation */

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
        double BendingAngle=atGetOptionalDouble(ElemData,"BendingAngle",0.0); check_error();
        double EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        double ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        double Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        int FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        int FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        int FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error();
        int FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        double gK=atGetOptionalDouble(ElemData,"gK", 0.0); check_error();
        double *R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        double *R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        double *T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        double *T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        double *EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        double *RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        double *KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();

        if (NumIntSteps == 0) {
            for (int i=MaxOrder; i>=0; i--) {
                if ((PolynomA[i] != 0.0) || (PolynomB[i] != 0.0)) {
                    atError("NumIntSteps == 0 not allowed with multipoles"); check_error();
                }
            }
        }

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        ExactSectorBend(r_in, Length, BendingAngle, PolynomA, PolynomB,
            MaxOrder, NumIntSteps, EntranceAngle, ExitAngle,
            FringeBendEntrance, FringeBendExit,
            FringeQuadEntrance, FringeQuadExit,
            gK,
            T1, T2, R1, R2, RApertures, EApertures,
            KickAngle, Scaling, num_particles);
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
            plhs[1] = mxCreateCellMatrix(13, 1);
            mxSetCell(plhs[1], i1++, mxCreateString("FringeBendEntrance"));
            mxSetCell(plhs[1], i1++, mxCreateString("FringeBendExit"));
            mxSetCell(plhs[1], i1++, mxCreateString("FringeQuadEntrance"));
            mxSetCell(plhs[1], i1++, mxCreateString("FringeQuadExit"));
            mxSetCell(plhs[1], i1++, mxCreateString("gK"));
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
