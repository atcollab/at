#include "mex.h"
#include "elempass.h"

void DriftPass(double *r_in, double le, int num_particles)
/* le - physical length
   r_in - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{
    double *r6;
    int c;
	double p_norm, NormL;
    
	for (c = 0; c<num_particles; c++) { /*Loop over particles  */
        r6 = r_in+c*6;
        if(!mxIsNaN(r6[0])) {
            p_norm = 1/(1+r6[4]);
            NormL = le*p_norm;
            r6[0] += NormL*r6[1];
            r6[2] += NormL*r6[3];
            r6[5] += NormL*p_norm*(r6[1]*r6[1]+r6[3]*r6[3])/2;
        }
    }
}

#ifndef NOMEX

#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
                             double *r_in, int num_particles, int mode)
#define NUM_FIELDS_2_REMEMBER 1
{
	double le;
    
	switch(mode) {
	    case MAKE_LOCAL_COPY: 	/* Find field numbers first
                                 Save a list of field number in an array
                                 and make returnptr point to that array
                                 */
            FieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            
            /*  Populate */
            
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "Length");
            /* Fall through next section... */
            
	    case USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
                                 The second argument ponter to the array of field
                                 numbers is previously created with
                                 QuadLinPass( ..., MAKE_LOCAL_COPY)
                                 */
            le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
            break;
	    default:
            mexErrMsgTxt("No match for calling mode in function QuadMPoleFringePass\n");
	}
	DriftPass(r_in, le, num_particles);
	return FieldNumbers;
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        double le = mxGetScalar(GetRequiredField(prhs[0], "Length"));
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        DriftPass(r_in,le, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(1,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif
