#include "at.h"
#include "atlalib.c"


#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656



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
    bool useT1, useT2, useR1, useR2;
    double SL, L1, L2, K1, K2;
    SL = le/num_int_steps;
    L1 = SL*DRIFT1;
    L2 = SL*DRIFT2;
    K1 = SL*KICK1;
    K2 = SL*KICK2;
    
    if(T1==NULL)
        useT1=false;
    else
        useT1=true;
    
    if(T2==NULL)
        useT2=false;
    else
        useT2=true;
    
    if(R1==NULL)
        useR1=false;
    else
        useR1=true;
    
    if(R2==NULL)
        useR2=false;
    else
        useR2=true;
    for(c = 0;c<num_particles;c++)	/*Loop over particles  */
    {	r6 = r+c*6;
        if(!mxIsNaN(r6[0]))
        {
            /*  misalignment at entrance  */
            if(useT1)
                ATaddvv(r6,T1);
            if(useR1)
                ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /*  integrator  */
            for(m=0; m < num_int_steps; m++) /*  Loop over slices */
            {	r6 = r+c*6;
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
            if(useR2)
                ATmultmv(r6,R2);
            if(useT2)
                ATaddvv(r6,T2);
        }
    }
}

#ifdef PYAT

#include "pyutils.c"

int atpyPass(double *rin, int num_particles, PyObject *element, struct parameters *param)
{
    PyErr_Clear();
    double *t1 = numpy_get_double_array(element, "T1");     /* Optional arguments */
    double *t2 = numpy_get_double_array(element, "T2");
    double *r1 = numpy_get_double_array(element, "R1");
    double *r2 = numpy_get_double_array(element, "R2");
    double *RApertures = numpy_get_double_array(element, "RApertures");
    double *EApertures = numpy_get_double_array(element, "EApertures");
    double length = py_get_double(element, "Length");       /* Mandatory arguments */
    long max_order = py_get_long(element, "MaxOrder");
    long num_int_steps = py_get_long(element, "NumIntSteps");
    double *polyA = numpy_get_double_array(element, "PolynomA");
    double *polyB = numpy_get_double_array(element, "PolynomB");
    if (PyErr_Occurred())
        return -1;
    else {
        StrMPoleSymplectic4Pass(rin, length, polyA, polyB, (int)max_order, (int)num_int_steps, t1, t2, r1, r2,
            RApertures, EApertures, num_particles);
        return 0;
    }
}

#endif /*PYAT*/

#ifdef MATLAB_MEX_FILE

#include "mex.h"
#include "elempass.h"

#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
        double *r_in, int num_particles, int mode)
#define NUM_FIELDS_2_REMEMBER 11
{	
    double *A , *B;
    double  *pr1, *pr2, *pt1, *pt2, *RApertures, *EApertures;
    int max_order, num_int_steps;
    double le;
    int *returnptr;
        
    switch(mode)
    {	case NO_LOCAL_COPY:	/* NOT used in AT1.3 Get fields by names from MATLAB workspace   */
        {
            
        }	break;
        
        case MAKE_LOCAL_COPY: 	/* Find field numbers first
         * Save a list of field number in an array
         * and make returnptr point to that array
         */
            /* Allocate memory for integer array of
             * field numbers for faster futurereference
             */
            FieldNumbers = (int *) mxCalloc(NUM_FIELDS_2_REMEMBER, sizeof(int));
            
            /*  Populate */
            
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "PolynomA");
            FieldNumbers[1] = GetRequiredFieldNumber(ElemData, "PolynomB");
            FieldNumbers[2] = GetRequiredFieldNumber(ElemData, "MaxOrder");
            FieldNumbers[3] = GetRequiredFieldNumber(ElemData, "NumIntSteps");
            FieldNumbers[4] = GetRequiredFieldNumber(ElemData, "Length");
            
            FieldNumbers[5] = mxGetFieldNumber(ElemData,"R1");
            FieldNumbers[6] = mxGetFieldNumber(ElemData,"R2");
            FieldNumbers[7] = mxGetFieldNumber(ElemData,"T1");
            FieldNumbers[8] = mxGetFieldNumber(ElemData,"T2");
            FieldNumbers[9] = mxGetFieldNumber(ElemData,"RApertures");
            FieldNumbers[10] = mxGetFieldNumber(ElemData,"EApertures");
            /* Fall through next section... */
            
        case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
         * The second argument ponter to the array of field
         * numbers is previously created with
         * StrMPoleSymplectic4Pass( ..., MAKE_LOCAL_COPY)
         *
         */
            A = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
            B = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
            max_order = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
            num_int_steps = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
            le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
            
            /* Optional fields */
            pr1 = (FieldNumbers[5] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[5])) : NULL;
            pr2 = (FieldNumbers[6] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[6])) : NULL;
            pt1 = (FieldNumbers[7] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[7])) : NULL;
            pt2 = (FieldNumbers[8] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[8])) : NULL;
            RApertures = (FieldNumbers[9] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[9])) : NULL;
            EApertures = (FieldNumbers[10] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[10])) : NULL;
            break;
        default:
            mexErrMsgTxt("No match for calling mode in function StrMPoleSymplectic4Pass\n");
            
    }
    
    StrMPoleSymplectic4Pass(r_in, le, A, B, max_order, num_int_steps,
            pt1, pt2, pr1, pr2, RApertures, EApertures, num_particles);
    return FieldNumbers;
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        double *pr1, *pr2, *pt1, *pt2, *RApertures, *EApertures;
        mxArray *tmpmxptr;
        
        double *A = mxGetPr(GetRequiredField(prhs[0], "PolynomA"));
        double *B = mxGetPr(GetRequiredField(prhs[0], "PolynomB"));
        int max_order = (int) mxGetScalar(GetRequiredField(prhs[0], "MaxOrder"));
        int num_int_steps = (int) mxGetScalar(GetRequiredField(prhs[0], "NumIntSteps"));
        double le = mxGetScalar(GetRequiredField(prhs[0], "Length"));
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        
        /* Optional arguments */
        tmpmxptr = mxGetField(prhs[0],0,"R1");
        pr1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"R2");
        pr2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T1");
        pt1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T2");
        pt2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"RApertures");
        RApertures = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"EApertures");
        EApertures = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        
        StrMPoleSymplectic4Pass(r_in, le, A, B, max_order, num_int_steps,
                pt1, pt2, pr1, pr2, RApertures, EApertures, num_particles);
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
