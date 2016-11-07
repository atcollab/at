/* DriftPass.c 
   Accelerator Toolbox 
   Revision 6/26/00
   A.Terebilo terebilo@ssrl.slac.stanford.edu
*/

#include "at.h"
#include "atlalib.c"

void Matrix66Pass(double *r, const double *M,
        const double *T1, const double *T2,
        const double *R1, const double *R2, int num_particles)

{
    double *r6;
    int c;
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        r6 = r+c*6;
        if (!mxIsNaN(r6[0])) {
            if (T1 != NULL) ATaddvv(r6, T1);
            if (R1 != NULL) ATmultmv(r6, R1);
            ATmultmv(r6, M);
            if (R2 != NULL) ATmultmv(r6, R2);
            if (T2 != NULL) ATaddvv(r6, T2);
        }
	}
}

MODULE_DEF(Matrix66Pass)        /* Dummy module initialisation */

#ifdef MATLAB_MEX_FILE

#include "elempass.h"
#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
        double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 5
{
    double *M;
    double  *pr1, *pr2, *pt1, *pt2;
    
    switch(mode) {
        case MAKE_LOCAL_COPY: 	/* Find field numbers first
         * Save a list of field number in an array
         * and make returnptr point to that array
         */
            FieldNumbers = (int *) mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "M66");
            
            
            FieldNumbers[1] = mxGetFieldNumber(ElemData,"R1");
            FieldNumbers[2] = mxGetFieldNumber(ElemData,"R2");
            FieldNumbers[3] = mxGetFieldNumber(ElemData,"T1");
            FieldNumbers[4] = mxGetFieldNumber(ElemData,"T2");
            /* Fall through next section... */
            
        case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
         * The second argument ponter to the array of field
         * numbers is previously created with
         * QuadLinPass( ..., MAKE_LOCAL_COPY)
         */
            M = mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[0]));
            
            /* Optional fields */
            pr1 = (FieldNumbers[1] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[1])) : NULL;
            pr2 = (FieldNumbers[2] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[2])) : NULL;
            pt1 = (FieldNumbers[3] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[3])) : NULL;
            pt2 = (FieldNumbers[4] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[4])) : NULL;
            break;
            
        default:
            mexErrMsgTxt("No match for calling mode in function QuadMPoleFringePass\n");
            
    }
    Matrix66Pass(r_in, M, pt1, pt2, pr1, pr2, num_particles);
    return FieldNumbers;
}


void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        double  *pr1, *pr2, *pt1, *pt2;
        mxArray *tmpmxptr;
        
        double *M = mxGetPr(GetRequiredField(prhs[0], "M66"));
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
 		
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        Matrix66Pass(r_in, M, pt1, pt2, pr1, pr2, num_particles);	
	}
    else if (nrhs == 0) {
        /* list of required fields */
	    plhs[0] = mxCreateCellMatrix(1,1);
	    mxSetCell(plhs[0],0,mxCreateString("M66"));
	    if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(4,1);
            mxSetCell(plhs[1],0,mxCreateString("T1"));
            mxSetCell(plhs[1],1,mxCreateString("T2"));
            mxSetCell(plhs[1],2,mxCreateString("R1"));
            mxSetCell(plhs[1],3,mxCreateString("R2"));
	    }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif
