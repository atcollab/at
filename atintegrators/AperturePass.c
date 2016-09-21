#include <math.h>
#include "at.h"

static void markaslost(double *r6,int idx)
{
    r6[idx] = mxGetInf();
}

void AperturePass(double *r_in, double *limits, int num_particles)
{
    /* Checks X and Y of each input 6-vector and marks the corresponding element in
     * lossflag array with 0 if X,Y are exceed the limits given by limitsptr array
     * limitsptr has 4 elements: (MinX, MaxX, MinY, MaxY)
     */
    int c;
    double *r6;
    for (c=0; c<num_particles; c++) {
        r6 = r_in+c*6;
        if (!mxIsNaN(r6[0])) { /*  check if this particle is already marked as lost */
            /* check limits for X position */
            if (r6[0]<limits[0] || r6[0]>limits[1])      markaslost(r6,0);
            else if (r6[2]<limits[2] || r6[2]>limits[3]) markaslost(r6,2);
        }
    }
}

#ifdef PYAT

#include "pyutils.c"

int atpyPass(double *rin, int num_particles, PyObject *element, struct parameters *param)
{
    double *limits = numpy_get_double_array(element, "Limits", false); /* Mandatory arguments */
    if (PyErr_Occurred()) {
        return -1;
    } else {
        AperturePass(rin, limits, num_particles);
        return 0;
    }
}

#endif /*PYAT*/

#ifdef MATLAB_MEX_FILE

#include "mex.h"
#include "elempass.h"
#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
        double *r_in, int num_particles, int mode)
        
#define NUM_FIELDS_2_REMEMBER 1
{
    double *limits;
    
    switch(mode) {
        case MAKE_LOCAL_COPY:
            /* Find field numbers first
             * Save a list of field number in an array
             * and make returnptr point to that array
             */
            FieldNumbers = (int *) mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "Limits");
            
            /* Fall through next section... */
            
        case USE_LOCAL_COPY:
            /* Get fields from MATLAB using field numbers
             * The second argument pointer to the array of field
             * numbers is previously created with
             * QuadLinPass( ..., MAKE_LOCAL_COPY)
             */
            limits = mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[0]));
            break;
            
        default:
            mexErrMsgTxt("No match for calling mode in function AperturePass\n");
    }
    
    AperturePass(r_in, limits, num_particles);
    
    return FieldNumbers;
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        
        double *limits = mxGetPr(GetRequiredField(prhs[0], "Limits"));
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        AperturePass(r_in,limits, num_particles);
    }
    else if (nrhs == 0) {
        /* return list of required fields */
        plhs[0] = mxCreateCellMatrix(1,1);
        mxSetCell(plhs[0],0,mxCreateString("Limits"));
        if(nlhs>1) {
            /* Required and optional fields */
            plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/
