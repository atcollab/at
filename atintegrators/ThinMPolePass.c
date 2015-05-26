/* ThinMPolePass
   Accelerator Toolbox 
   A.Terebilo terebilo@ssrl.slac.stanford.edu
*/

#include "mex.h"
#include "elempass.h"
#include "atlalib.c"
#include "driftkick.c"

void ThinMPolePass(double *r, double *A, double *B, int max_order,
        double bax, double bay,
        double *T1, double *T2,
        double *R1, double *R2, int num_particles)
{
    int c;
    double *r6;
    for (c = 0;c<num_particles;c++)	{ /* Loop over particles */
        r6 = r+c*6;
        if(!mxIsNaN(r6[0])) {
            /*  misalignment at entrance  */
            if (T1 != NULL) ATaddvv(r6,T1);
            if (R1 != NULL) ATmultmv(r6,R1);
            strthinkick(r6, A, B, 1.0, max_order);
            r6[1] += bax*r6[4];
            r6[3] -= bay*r6[4];
            r6[5] -= bax*r6[0]-bay*r6[2]; /* Path lenghtening */
            /*  misalignment at exit  */
            if (R2 != NULL) ATmultmv(r6,R2);
            if (T2 != NULL) ATaddvv(r6,T2);
        }
    }
}

#ifndef NOMEX

#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
								double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 8

{
    double *A , *B, bax, bay;
    double  *pr1, *pr2, *pt1, *pt2;
    int max_order;
    
    switch(mode) {
        /* case NO_LOCAL_COPY:	Not used in AT1.3 */
        
        case MAKE_LOCAL_COPY: 	/* Find field numbers first
                                    Save a list of field number in an array
                                    and make returnptr point to that array
                                */
            /* Allocate memory for integer array of
             * field numbers for faster futurereference
             */
            
            FieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            
            /* Populate */
            
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "PolynomA");
            FieldNumbers[1] = GetRequiredFieldNumber(ElemData, "PolynomB");
            FieldNumbers[2] = GetRequiredFieldNumber(ElemData, "MaxOrder");
            
            /* Optional fields */
            FieldNumbers[3] = mxGetFieldNumber(ElemData,"BendingAngle");
            FieldNumbers[4] = mxGetFieldNumber(ElemData,"R1");
            FieldNumbers[5] = mxGetFieldNumber(ElemData,"R2");
            FieldNumbers[6] = mxGetFieldNumber(ElemData,"T1");
            FieldNumbers[7] = mxGetFieldNumber(ElemData,"T2");
            
        case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
         * The second argument ponter to the array of field
         * numbers is previously created with
         * QuadLinPass( ..., MAKE_LOCAL_COPY)
         */
            
            A = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
            B = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
            max_order = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
            if (FieldNumbers[3] >= 0) {
                mxArray *tmpmxptr = mxGetFieldByNumber(ElemData,0,FieldNumbers[3]);
                double *tmpdblptr = mxGetPr(tmpmxptr);
                bax = tmpdblptr[0];
                bay = (mxGetNumberOfElements(tmpmxptr)>1) ? tmpdblptr[1] : 0.0;
            }
            else {
                bax = 0.0;
                bay = 0.0;
            }
            /* Optional fields */
            pr1 = (FieldNumbers[4] < 0) ? NULL : mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
            pr2 = (FieldNumbers[5] < 0) ? NULL : mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
            pt1 = (FieldNumbers[6] < 0) ? NULL : mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
            pt2 = (FieldNumbers[7] < 0) ? NULL : mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));
            break;
        default:
            mexErrMsgTxt("No match for calling mode in function ThinMPolePass\n");
    }
    
    ThinMPolePass(r_in,A, B, max_order, bax, bay, pt1, pt2, pr1, pr2, num_particles);
    
    return FieldNumbers;
    
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        double bax, bay, *pr1, *pr2, *pt1, *pt2 ;
        mxArray *tmpmxptr;
        
        double *A = mxGetPr(GetRequiredField(prhs[0], "PolynomA"));
        double *B = mxGetPr(GetRequiredField(prhs[0], "PolynomB"));
        int max_order = (int) mxGetScalar(GetRequiredField(prhs[0], "MaxOrder"));
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
                
        /* Optional arguments */
        tmpmxptr = mxGetField(prhs[0],0,"BendingAngle");
        if (tmpmxptr) {
            double *tmpdblptr = mxGetPr(tmpmxptr);
            bax = tmpdblptr[0];
            bay = (mxGetNumberOfElements(tmpmxptr)>1) ? tmpdblptr[1] : 0.0;
        }
        else {
            bax = 0;
            bay = 0;
        }
        tmpmxptr = mxGetField(prhs[0],0,"R1");
        pr1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"R2");
        pr2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T1");
        pt1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T2");
        pt2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;

        /* ALLOCATE memory for the output array of the same size as the input */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        ThinMPolePass(r_in, A, B, max_order, bax, bay,
                pt1, pt2, pr1, pr2, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(3,1);
        mxSetCell(plhs[0],0,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],1,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],2,mxCreateString("MaxOrder"));
        if(nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(5,1);
            mxSetCell(plhs[1],0,mxCreateString("BendingAngle"));
            mxSetCell(plhs[1],1,mxCreateString("T1"));
            mxSetCell(plhs[1],2,mxCreateString("T2"));
            mxSetCell(plhs[1],3,mxCreateString("R1"));
            mxSetCell(plhs[1],4,mxCreateString("R2"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
    
}
#endif
