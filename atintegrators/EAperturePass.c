#include "mex.h"
#include <math.h>
#include "elempass.h"

static void markaslost(double *r6,int idx)
{
    r6[idx] = mxGetInf();
}

void EAperturePass(double *r_in, double *axesptr, int num_particles)
{
    /*  Checks X and Y of each input 6-vector and marks the corresponding element in
     * lossflag array with 0 if X,Y are exceed the limits given by axesptr array
     */
    int c;
    double *r6;
    for (c=0;c<num_particles;c++) {
        r6 = r_in+c*6;
        if (!mxIsNaN(r6[0])) { /*  check if this particle is already marked as lost */
            register double xnorm = r6[0]/axesptr[0];
            register double znorm = r6[2]/axesptr[1];
            if ((xnorm*xnorm + znorm*znorm) >= 1) markaslost(r6,0);
        }
    }
}

#ifndef NOMEX

#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
        double *r_in, int num_particles, int mode)
        
#define NUM_FIELDS_2_REMEMBER 1
{
    double *axes;
    
    switch(mode) {
        case MAKE_LOCAL_COPY:
            /* Find field numbers first
             * Save a list of field number in an array
             * and make returnptr point to that array
             */
            FieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "Axes");
            
            /* Fall through next section... */
        case	USE_LOCAL_COPY:
            /* Get fields from MATLAB using field numbers
             * The second argument ponter to the array of field
             * numbers is previously created with
             * QuadLinPass( ..., MAKE_LOCAL_COPY)
             */
            axes = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
            break;
            
        default:
            mexErrMsgTxt("No match for calling mode in function EAperturePass\n");
    }
    
    EAperturePass(r_in, axes, num_particles);
    
    return FieldNumbers;
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        
        double *axes = mxGetPr(GetRequiredField(prhs[0], "Axes"));
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        EAperturePass(r_in, axes,num_particles);
    }
    else if (nrhs == 0) {
        /* return list of required fields */
        plhs[0] = mxCreateCellMatrix(1,1);
        mxSetCell(plhs[0],0,mxCreateString("Axes"));
        if(nlhs>1) {
            /* Required and optional fields */
            plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*NOMEX*/
