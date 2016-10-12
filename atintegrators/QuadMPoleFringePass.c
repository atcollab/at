
#include "atelem.c"
#include "atlalib.c"
#include "driftkick.c"		/* fastdrift.c, strthinkick.c */
#include "quadfringe.c"		/* QuadFringePassP, QuadFringePassN */

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

struct elem
{
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    /* Optional fields */
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *fringeIntM0;
    double *fringeIntP0;
};

void QuadMPoleFringePass(double *r, double le, const double *A, const double *B,
        int max_order, int num_int_steps,
        const double *T1, const double *T2,
        const double *R1, const double *R2,
        double *fringeIntM0,  /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
        double *fringeIntP0,  /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */
        double *RApertures, double *EApertures, int num_particles)
{
    double *r6;
    int c, m;
    double p_norm, NormL1, NormL2;
    bool useLinFrEle = (fringeIntM0 != NULL && fringeIntP0 != NULL);
    double SL = le/num_int_steps;
    double L1 = SL*DRIFT1;
    double L2 = SL*DRIFT2;
    double K1 = SL*KICK1;
    double K2 = SL*KICK2;
    double *fringeIntM, *fringeIntP;
    double delta;
    double inFringe; 	/*  for linear fringe field, from elegant.
     * The argument inFringe is a flag:
     * -1 for going into the magnet and +1 for going out. */
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        r6 = r+c*6;
        if (!atIsNaN(r6[0])) {
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6, T1);
            if (R1) ATmultmv(r6, R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            if (useLinFrEle) /*Linear fringe fields from elegant*/
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
                /* nonlinear fringe field */
                QuadFringePassP(r6,B[1]);
            
            /*  integrator  */
            p_norm = 1/(1+r6[4]);
            NormL1 = L1*p_norm;
            NormL2 = L2*p_norm;
            for (m=0; m < num_int_steps; m++) { /*  Loop over slices */
                fastdrift(r6, NormL1);
                strthinkick(r6, A, B,  K1, max_order);
                fastdrift(r6, NormL2);
                strthinkick(r6, A, B, K2, max_order);
                fastdrift(r6, NormL2);
                strthinkick(r6, A, B,  K1, max_order);
                fastdrift(r6, NormL1);
            }
            if (useLinFrEle) /*Linear fringe fields from elegant*/
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
                /* nonlinear fringe field */
                QuadFringePassN(r6,B[1]);
            
            /* Check physical apertures at the exit of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* Misalignment at exit */
            if (R2) ATmultmv(r6, R2);
            if (T2) ATaddvv(r6, T2);
        }
    }
}



#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Length;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        Elem->NumIntSteps=NumIntSteps;
        /*optional fields*/
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->fringeIntM0=fringeIntM0;
        Elem->fringeIntP0=fringeIntP0;
    }
    QuadMPoleFringePass(r_in,Elem->Length,Elem->PolynomA,Elem->PolynomB,
            Elem->MaxOrder,Elem->NumIntSteps,Elem->T1,Elem->T2,Elem->R1,Elem->R2,
            Elem->fringeIntM0,Elem->fringeIntP0,Elem->RApertures,Elem->EApertures,num_particles);
    return Elem;
}

void initQuadMPoleFringePass(void) {};
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2 ) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        QuadMPoleFringePass(r_in,Length,PolynomA,PolynomB,MaxOrder,NumIntSteps,
                T1,T2,R1,R2,fringeIntM0,fringeIntP0,RApertures,EApertures,num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(5,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],2,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],3,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],4,mxCreateString("NumIntSteps"));

        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(8,1);
            mxSetCell(plhs[1],0,mxCreateString("T1"));
            mxSetCell(plhs[1],1,mxCreateString("T2"));
            mxSetCell(plhs[1],2,mxCreateString("R1"));
            mxSetCell(plhs[1],3,mxCreateString("R2"));
            mxSetCell(plhs[1],4,mxCreateString("RApertures"));
            mxSetCell(plhs[1],5,mxCreateString("EApertures"));
            mxSetCell(plhs[1],6,mxCreateString("fringeIntM0"));
            mxSetCell(plhs[1],7,mxCreateString("fringeIntP0"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
