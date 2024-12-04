#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "atphyslib.c"
#include "driftkick.c"		/* fastdrift and bndthinkick */
#include "quadfringe.c"		/* QuadFringePassP, QuadFringePassN */

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
    int FringeBendEntrance;
    int FringeBendExit;
    double FringeInt1;
    double FringeInt2;
    double FullGap;
    double Scaling;
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
};

void BndMPoleSymplectic4Pass(double *r, double le, double irho, double *A, double *B,
        int max_order, int num_int_steps,
        double entrance_angle, double exit_angle,
        int FringeBendEntrance, int FringeBendExit,
        double fint1, double fint2, double gap,
        int FringeQuadEntrance, int FringeQuadExit,
        double *fringeIntM0,  /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
        double *fringeIntP0,  /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double *KickAngle, double scaling, int num_particles)
{
    double SL = le/num_int_steps;
    double L1 = SL*DRIFT1;
    double L2 = SL*DRIFT2;
    double K1 = SL*KICK1;
    double K2 = SL*KICK2;
    bool useLinFrEleEntrance = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadEntrance==2);
    bool useLinFrEleExit = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadExit==2);
    double B0 = B[0];
    double A0 = A[0];

    if (KickAngle) {   /* Convert corrector component to polynomial coefficients */
        B[0] -= sin(KickAngle[0])/le;
        A[0] += sin(KickAngle[1])/le;
    }

    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
    shared(r,num_particles,R1,T1,R2,T2,RApertures,EApertures,\
    irho,gap,A,B,L1,L2,K1,K2,max_order,num_int_steps,scaling,\
    FringeBendEntrance,entrance_angle,fint1,FringeBendExit,exit_angle,fint2,\
    FringeQuadEntrance,useLinFrEleEntrance,FringeQuadExit,useLinFrEleExit,fringeIntM0,fringeIntP0)
    for (int c = 0; c<num_particles; c++) { /* Loop over particles */
        double *r6 = r + 6*c;
        if (!atIsNaN(r6[0])) {
            int m;
            double p_norm, NormL1, NormL2;
            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, scaling);
            p_norm = 1.0/(1.0+r6[4]);
            NormL1 = L1*p_norm;
            NormL2 = L2*p_norm;
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* edge focus */
            edge_fringe_entrance(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance);
            /* quadrupole gradient fringe entrance*/
            if (FringeQuadEntrance && B[1]!=0) {
                if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
                    linearQuadFringeElegantEntrance(r6, B[1], fringeIntM0, fringeIntP0);
                else
                    QuadFringePassP(r6, B[1]);
            }
            /* integrator */
            for (m=0; m < num_int_steps; m++) { /* Loop over slices */
                fastdrift(r6, NormL1);
                bndthinkick(r6, A, B, K1, irho, max_order);
                fastdrift(r6, NormL2);
                bndthinkick(r6, A, B, K2, irho, max_order);
                fastdrift(r6, NormL2);
                bndthinkick(r6, A, B, K1, irho, max_order);
                fastdrift(r6, NormL1);
            }
            /* quadrupole gradient fringe */
            if (FringeQuadExit && B[1]!=0) {
                if (useLinFrEleExit) /*Linear fringe fields from elegant*/
                    linearQuadFringeElegantExit(r6, B[1], fringeIntM0, fringeIntP0);
                else
                    QuadFringePassN(r6, B[1]);
            }
            /* edge focus */
            edge_fringe_exit(r6, irho, exit_angle, fint2, gap, FringeBendExit);
            /* Check physical apertures at the exit of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);
            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, 1.0/scaling);
        }
    }
    if (KickAngle) {  /* Remove corrector component in polynomial coefficients */
        B[0] = B0;
        A[0] = A0;
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double irho;
    if (!Elem) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, Scaling,
                FringeInt1, FringeInt2;
        int MaxOrder, NumIntSteps,  FringeBendEntrance, FringeBendExit,
                FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0, *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error();
        FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();

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
        Elem->FringeBendEntrance=FringeBendEntrance;
        Elem->FringeBendExit=FringeBendExit;
        Elem->FullGap=FullGap;
        Elem->Scaling=Scaling;
        Elem->FringeInt1=FringeInt1;
        Elem->FringeInt2=FringeInt2;
        Elem->FringeQuadEntrance=FringeQuadEntrance;
        Elem->FringeQuadExit=FringeQuadExit;
        Elem->fringeIntM0=fringeIntM0;
        Elem->fringeIntP0=fringeIntP0;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->KickAngle=KickAngle;
    }
    irho = Elem->BendingAngle/Elem->Length;
    BndMPoleSymplectic4Pass(r_in, Elem->Length, irho, Elem->PolynomA, Elem->PolynomB,
            Elem->MaxOrder, Elem->NumIntSteps, Elem->EntranceAngle, Elem->ExitAngle,
            Elem->FringeBendEntrance,Elem->FringeBendExit,
            Elem->FringeInt1, Elem->FringeInt2, Elem->FullGap,
            Elem->FringeQuadEntrance, Elem->FringeQuadExit,
            Elem->fringeIntM0, Elem->fringeIntP0,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            Elem->RApertures, Elem->EApertures,
            Elem->KickAngle, Elem->Scaling, num_particles);
    return Elem;
}

MODULE_DEF(BndMPoleSymplectic4Pass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, Scaling,
                FringeInt1, FringeInt2;
        int MaxOrder, NumIntSteps, FringeBendEntrance, FringeBendExit,
                FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0, *KickAngle;
        double irho;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgTxt("Second argument must be a 6 x N matrix");

        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error();
        FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();
        irho = BendingAngle/Length;

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        BndMPoleSymplectic4Pass(r_in, Length, irho, PolynomA, PolynomB,
            MaxOrder, NumIntSteps, EntranceAngle, ExitAngle,
            FringeBendEntrance, FringeBendExit,
            FringeInt1, FringeInt2, FullGap,
            FringeQuadEntrance, FringeQuadExit,
            fringeIntM0, fringeIntP0,
            T1, T2, R1, R2, RApertures, EApertures,
            KickAngle, Scaling, num_particles);
    } else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(8,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("BendingAngle"));
        mxSetCell(plhs[0],2,mxCreateString("EntranceAngle"));
        mxSetCell(plhs[0],3,mxCreateString("ExitAngle"));
        mxSetCell(plhs[0],4,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],5,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],6,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],7,mxCreateString("NumIntSteps"));

        if (nlhs>1) {    /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(17,1);
            mxSetCell(plhs[1],0,mxCreateString("FullGap"));
            mxSetCell(plhs[1],1,mxCreateString("FringeInt1"));
            mxSetCell(plhs[1],2,mxCreateString("FringeInt2"));
            mxSetCell(plhs[1],3,mxCreateString("FringeBendEntrance"));
            mxSetCell(plhs[1],4,mxCreateString("FringeBendExit"));
            mxSetCell(plhs[1],5,mxCreateString("FringeQuadEntrance"));
            mxSetCell(plhs[1],6,mxCreateString("FringeQuadExit"));
            mxSetCell(plhs[1],7,mxCreateString("fringeIntM0"));
            mxSetCell(plhs[1],8,mxCreateString("fringeIntP0"));
            mxSetCell(plhs[1],9,mxCreateString("T1"));
            mxSetCell(plhs[1],10,mxCreateString("T2"));
            mxSetCell(plhs[1],11,mxCreateString("R1"));
            mxSetCell(plhs[1],12,mxCreateString("R2"));
            mxSetCell(plhs[1],13,mxCreateString("RApertures"));
            mxSetCell(plhs[1],14,mxCreateString("EApertures"));
            mxSetCell(plhs[1],15,mxCreateString("KickAngle"));
            mxSetCell(plhs[1],16,mxCreateString("FieldScaling"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
