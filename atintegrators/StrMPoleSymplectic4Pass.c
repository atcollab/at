
#include "atelem.c"
#include "atlalib.c"

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
};

static void strthinkick(double* r, double* A, double* B, double L, int max_order)
/*****************************************************************************
 * Calculate and apply a multipole kick to a 6-dimentional
 * phase space vector in a straight element ( quadrupole)
 *
 * IMPORTANT !!!
 * The reference coordinate system is straight but the field expansion may still
 * contain dipole terms: PolynomA(1), PolynomB(1) - in MATLAB notation,
 * A[0], B[0] - C,C++ notation
 *
 *
 * Note: in the US convention the transverse multipole field is written as:
 *
 * max_order+1
 * ----
 * \                       n-1
 * (B + iB  )/ B rho  =  >   (ia  + b ) (x + iy)
 * y    x            /       n    n
 * ----
 * n=1
 * is a polynomial in (x,y) with the highest order = MaxOrder
 *
 *
 * Using different index notation
 *
 * max_order
 * ----
 * \                       n
 * (B + iB  )/ B rho  =  >   (iA  + B ) (x + iy)
 * y    x            /       n    n
 * ----
 * n=0
 *
 * A,B: i=0 ... max_order
 * [0] - dipole, [1] - quadrupole, [2] - sextupole ...
 * units for A,B[i] = 1/[m]^(i+1)
 * Coeficients are stroed in the PolynomA, PolynomB field of the element
 * structure in MATLAB
 *
 * A[i] (C++,C) =  PolynomA(i+1) (MATLAB)
 * B[i] (C++,C) =  PolynomB(i+1) (MATLAB)
 * i = 0 .. MaxOrder
 *
 ******************************************************************************/
{  int i;
   double ReSum = B[max_order];
   double ImSum = A[max_order];
   double ReSumTemp;
   for(i=max_order-1;i>=0;i--)
   {   ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
       ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
       ReSum = ReSumTemp;
   }
   
   r[1] -=  L*ReSum;
   r[3] +=  L*ImSum;
}

void fastdrift(double* r, double NormL)

/*   NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
 * in the loop if momentum deviation (delta) does not change
 * such as in 4-th order symplectic integrator w/o radiation
 */

{   double dx = NormL*r[1];
    double dy = NormL*r[3];
    r[0]+= dx;
    r[2]+= dy;
    r[5]+= NormL*(r[1]*r[1]+r[3]*r[3])/(2*(1+r[4]));
}

void StrMPoleSymplectic4Pass(double *r, double le, double *A, double *B,
        int max_order, int num_int_steps,
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures, int num_particles)
{	int c,m;
    double norm, NormL1, NormL2;
    double *r6;
    double SL, L1, L2, K1, K2;
    SL = le/num_int_steps;
    L1 = SL*DRIFT1;
    L2 = SL*DRIFT2;
    K1 = SL*KICK1;
    K2 = SL*KICK2;
    
    for (c = 0;c<num_particles;c++)	{   /*Loop over particles  */
        r6 = r+c*6;
        if(!atIsNaN(r6[0])) {
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /*  integrator  */
            for (m=0; m < num_int_steps; m++) {  /*  Loop over slices */
             	r6 = r+c*6;
                norm = 1/(1+r6[4]);
                NormL1 = L1*norm;
                NormL2 = L2*norm;
                fastdrift(r6, NormL1);
                strthinkick(r6, A, B,  K1, max_order);
                fastdrift(r6, NormL2);
                strthinkick(r6, A, B, K2, max_order);
                fastdrift(r6, NormL2);
                strthinkick(r6, A, B,  K1, max_order);
                fastdrift(r6, NormL1);
            }
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
        double Length;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures;
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
    }
    StrMPoleSymplectic4Pass(r_in,Elem->Length,Elem->PolynomA,Elem->PolynomB,
            Elem->MaxOrder,Elem->NumIntSteps,Elem->T1,Elem->T2,Elem->R1,Elem->R2,
            Elem->RApertures,Elem->EApertures,num_particles);
    return Elem;
}

MODULE_DEF(StrMPoleSymplectic4Pass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures;
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
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        StrMPoleSymplectic4Pass(r_in,Length,PolynomA,PolynomB,MaxOrder,NumIntSteps,
                T1,T2,R1,R2,RApertures,EApertures,num_particles);
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
            plhs[1] = mxCreateCellMatrix(6,1);
            mxSetCell(plhs[1],0,mxCreateString("T1"));
            mxSetCell(plhs[1],1,mxCreateString("T2"));
            mxSetCell(plhs[1],2,mxCreateString("R1"));
            mxSetCell(plhs[1],3,mxCreateString("R2"));
            mxSetCell(plhs[1],4,mxCreateString("RApertures"));
            mxSetCell(plhs[1],5,mxCreateString("EApertures"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/
