/* IdentityPass.c 
   Accelerator Toolbox
   Revision 7/16/03
   A.Terebilo terebilo@slac.stanford.edu
*/

#include "at.h"
#include "atlalib.c"


void IdentityPass(double *r_in,
        const double *T1, const double *T2,
        const double *R1, const double *R2,
        const double *limits, const double *axesptr,
        int num_particles)
{	
    double *r6;
    int c;
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        r6 = r_in+c*6;
        if (!mxIsNaN(r6[0])) {
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6, T1);
            if (R1) ATmultmv(r6, R1);
			/* Check physical apertures */
			if (limits) checkiflostRectangularAp(r6,limits);
			if (axesptr) checkiflostEllipticalAp(r6,axesptr);
            /* Misalignment at exit */
            if (R2) ATmultmv(r6, R2);
            if (T2) ATaddvv(r6, T2);
        }
    }
}

#ifdef PYAT

#include "pyutils.c"

int atpyPass(double *rin, int num_particles, PyObject *element, struct parameters *param)
{
    double *t1 = numpy_get_double_array(element, "T1", true);     /* Optional arguments */
    double *t2 = numpy_get_double_array(element, "T2", true);
    double *r1 = numpy_get_double_array(element, "R1", true);
    double *r2 = numpy_get_double_array(element, "R2", true);
    double *RApertures = numpy_get_double_array(element, "RApertures", true);
    double *EApertures = numpy_get_double_array(element, "EApertures", true);
    if (PyErr_Occurred()) {
        return -1;
    } else {
        IdentityPass(rin, t1, t2, r1, r2, RApertures, EApertures, num_particles);
        return 0;
    }
}

#endif /*PYAT*/

#ifdef MATLAB_MEX_FILE

#include "elempass.h"
#include "mxutils.c"
#define NUM_FIELDS_2_REMEMBER 6

ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
                             double *r_in, int num_particles, int mode)
{
    double  *pr1, *pr2, *pt1, *pt2;
    double *limits, *axesptr;

    switch(mode) {
        case MAKE_LOCAL_COPY: 	/* Find field numbers first
                                 Save a list of field number in an array
                                 and make returnptr point to that array
                                 */
            /* Allocate memory for integer array of
             field numbers for faster futurereference
             */
            FieldNumbers = (int *) mxCalloc(NUM_FIELDS_2_REMEMBER, sizeof(int));
            
            /*  Populate */
            
            FieldNumbers[0] = mxGetFieldNumber(ElemData,"R1");
            FieldNumbers[1] = mxGetFieldNumber(ElemData,"R2");
            FieldNumbers[2] = mxGetFieldNumber(ElemData,"T1");
            FieldNumbers[3] = mxGetFieldNumber(ElemData,"T2");
			FieldNumbers[4] = mxGetFieldNumber(ElemData,"RApertures");
            FieldNumbers[5] = mxGetFieldNumber(ElemData,"EApertures");
            /* Fall through next section... */
            
        case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
                                 The second argument pointer to the array of field
                                 numbers is previously created with
                                 passFunction(..., MAKE_LOCAL_COPY)
                                 */
            
            /* Optional fields */
            pr1 = (FieldNumbers[0] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[0])) : NULL;
            pr2 = (FieldNumbers[1] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[1])) : NULL;
            pt1 = (FieldNumbers[2] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[2])) : NULL;
            pt2 = (FieldNumbers[3] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[3])) : NULL;
            limits = (FieldNumbers[4] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[4])) : NULL;
			axesptr = (FieldNumbers[5] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[5])) : NULL;
            
            break;
            
        default:
            mexErrMsgTxt("No match for calling mode in function IdentityPass\n");
    }
    
    IdentityPass(r_in, pt1, pt2, pr1, pr2, limits, axesptr, num_particles);
    return FieldNumbers;
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        double *pr1, *pr2, *pt1, *pt2;
        double *limits, *axesptr;
        mxArray *tmpmxptr;
        
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
        limits = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"EApertures");
        axesptr = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        IdentityPass(r_in, pt1, pt2, pr1, pr2, limits, axesptr, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(0,0);
        
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
#endif
