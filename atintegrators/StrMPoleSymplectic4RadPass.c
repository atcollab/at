#include "atelem.c"
#include "atlalib.c"
#include "driftkickrad.c"	/* drift6.c, strthinkickrad.c */
#include "quadfringe.c"		/* QuadFringePassP, QuadFringePassN */

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

#define SQR(X) ((X)*(X))

struct elem
{
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    double Energy;
    /* Optional fields */
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
};

void StrMPoleSymplectic4RadPass(double *r, double le, double *A, double *B,
        int max_order, int num_int_steps,
        int FringeQuadEntrance, int FringeQuadExit, /* 0 (no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) */
        double *fringeIntM0,  /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
        double *fringeIntP0,  /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double E0,
        int num_particles)
{	int c,m;
    double *r6;
    double SL, L1, L2, K1, K2;
    bool useLinFrEleEntrance = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadEntrance==2);
    bool useLinFrEleExit = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadExit==2);
    double *fringeIntM, *fringeIntP;  /*  for linear fringe field, from elegant. */
    double delta, inFringe; 	/*  for linear fringe field, from elegant. */
    SL = le/num_int_steps;
    L1 = SL*DRIFT1;
    L2 = SL*DRIFT2;
    K1 = SL*KICK1;
    K2 = SL*KICK2;
    
    for (c = 0;c<num_particles;c++)	{   /* Loop over particles  */
        r6 = r+c*6;
        if(!atIsNaN(r6[0])) {
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            if (FringeQuadEntrance && B[1]!=0)
                if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
                {
                    double R[6][6];
                    /* quadrupole linear fringe field, from elegant code */
                    inFringe=-1.0;
                    fringeIntM = fringeIntP0;
                    fringeIntP = fringeIntM0;
                    delta = r6[4];
                    /* determine first linear matrix for this delta */
                    quadPartialFringeMatrix(R, B[1]/(1+delta), inFringe, fringeIntM, 1);
                    r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
                    r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
                    r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
                    r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
                    /* nonlinear fringe field */
                    QuadFringePassP(r6,B[1]);   /*This is original AT code*/
                    /*Linear fringe fields from elegant*/
                    inFringe=-1.0;
                    /* determine and apply second linear matrix, from elegant code */
                    quadPartialFringeMatrix(R, B[1]/(1+delta), inFringe, fringeIntP, 2);
                    r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
                    r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
                    r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
                    r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
                }	/* end of elegant code*/
                else
                    QuadFringePassP(r6,B[1]);
            /* integrator */
            for (m=0; m < num_int_steps; m++) { /* Loop over slices */
             		r6 = r+c*6;
                    ATdrift6(r6,L1);
                    strthinkickrad(r6, A, B, K1, E0, max_order);
                    ATdrift6(r6,L2);
                    strthinkickrad(r6, A, B, K2, E0, max_order);
                    ATdrift6(r6,L2);
                    strthinkickrad(r6, A, B,  K1, E0, max_order);
                    ATdrift6(r6,L1);
            }
            if (FringeQuadExit && B[1]!=0)
                if (useLinFrEleExit) /*Linear fringe fields from elegant*/
                {
                double R[6][6];
                /* quadrupole linear fringe field, from elegant code */
                inFringe=1.0;
                fringeIntM = fringeIntM0;
                fringeIntP = fringeIntP0;
                delta = r6[4];
                /* determine first linear matrix for this delta */
                quadPartialFringeMatrix(R, B[1]/(1+delta), inFringe, fringeIntM, 1);
                r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
                r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
                r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
                r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
                /* nonlinear fringe field */
                QuadFringePassN(r6,B[1]);   /*This is original AT code*/
                /*Linear fringe fields from elegant*/
                inFringe=1.0;
                /* determine and apply second linear matrix, from elegant code */
                quadPartialFringeMatrix(R, B[1]/(1+delta), inFringe, fringeIntP, 2);
                r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
                r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
                r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
                r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
                }	/* end of elegant code*/
                else
                    QuadFringePassN(r6,B[1]);
            /* Check physical apertures at the exit of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2); 
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Length, Energy;
        int MaxOrder, NumIntSteps, FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
        FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0);
        FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0);
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        Elem->NumIntSteps=NumIntSteps;
        Elem->Energy=Energy;
        /*optional fields*/
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
    }
    StrMPoleSymplectic4RadPass(r_in,Elem->Length,Elem->PolynomA,Elem->PolynomB,
            Elem->MaxOrder,Elem->NumIntSteps,Elem->FringeQuadEntrance,
            Elem->FringeQuadExit,Elem->fringeIntM0,Elem->fringeIntP0,
            Elem->T1,Elem->T2,Elem->R1,Elem->R2,
            Elem->RApertures,Elem->EApertures,Elem->Energy,num_particles);
    return Elem;
}

MODULE_DEF(StrMPoleSymplectic4RadPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length, Energy;
        int MaxOrder, NumIntSteps, FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
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
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        StrMPoleSymplectic4RadPass(r_in,Length,PolynomA,PolynomB,MaxOrder,NumIntSteps,
                FringeQuadEntrance,FringeQuadExit,fringeIntM0,fringeIntP0,
                T1,T2,R1,R2,RApertures,EApertures,Energy,num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],2,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],3,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],4,mxCreateString("NumIntSteps"));
        mxSetCell(plhs[0],5,mxCreateString("Energy"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(10,1);
            mxSetCell(plhs[1],0,mxCreateString("FringeQuadEntrance"));
            mxSetCell(plhs[1],1,mxCreateString("FringeQuadExit"));  
            mxSetCell(plhs[1],2,mxCreateString("fringeIntM0"));
            mxSetCell(plhs[1],3,mxCreateString("fringeIntP0"));
            mxSetCell(plhs[1],4,mxCreateString("T1"));
            mxSetCell(plhs[1],5,mxCreateString("T2"));
            mxSetCell(plhs[1],6,mxCreateString("R1"));
            mxSetCell(plhs[1],7,mxCreateString("R2"));
            mxSetCell(plhs[1],8,mxCreateString("RApertures"));
            mxSetCell(plhs[1],9,mxCreateString("EApertures"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/

