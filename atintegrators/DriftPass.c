#include "at.h"
#include "atlalib.c"

void DriftPass(double *r_in, double le,
	       const double *T1, const double *T2,
	       const double *R1, const double *R2,
	       double *RApertures, double *EApertures,
	       int num_particles)
/* le - physical length
   r_in - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{
  double *r6;
  int c;

  for (c = 0; c<num_particles; c++) { /*Loop over particles  */
    r6 = r_in+c*6;
    if(!mxIsNaN(r6[0])) {
      /*  misalignment at entrance  */
      if (T1) ATaddvv(r6, T1);
      if (R1) ATmultmv(r6, R1);
      /* Check physical apertures at the entrance of the magnet */
      if (RApertures) checkiflostRectangularAp(r6,RApertures);
      if (EApertures) checkiflostEllipticalAp(r6,EApertures);
      ATdrift6(r6, le);
      /* Check physical apertures at the exit of the magnet */
      if (RApertures) checkiflostRectangularAp(r6,RApertures);
      if (EApertures) checkiflostEllipticalAp(r6,EApertures);
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
    PyErr_Clear();
    double *t1 = numpy_get_double_array(element, "T1");     /* Optional arguments */
    double *t2 = numpy_get_double_array(element, "T2");
    double *r1 = numpy_get_double_array(element, "R1");
    double *r2 = numpy_get_double_array(element, "R2");
    double *RApertures = numpy_get_double_array(element, "RApertures");
    double *EApertures = numpy_get_double_array(element, "EApertures");
    double length = py_get_double(element, "Length");       /* Mandatory arguments */
    if (PyErr_Occurred())
        return -1;
    else {
        DriftPass(rin, length, t1, t2, r1, r2, RApertures, EApertures, num_particles);
        return 0;
    }
}

#endif /*PYAT*/


#ifdef MATLAB_MEX_FILE

#include <mex.h>
#include "elempass.h"
#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
                             double *r_in, int num_particles, int mode)
#define NUM_FIELDS_2_REMEMBER 7
{
	double le;
    double  *pr1, *pr2, *pt1, *pt2;
    double *RApertures, *EApertures;

	switch(mode) {
	    case MAKE_LOCAL_COPY: 	/* Find field numbers first
                                 Save a list of field number in an array
                                 and make returnptr point to that array
                                 */
            FieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            
            /*  Populate */
            
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "Length");

            FieldNumbers[1] = mxGetFieldNumber(ElemData,"R1");
            FieldNumbers[2] = mxGetFieldNumber(ElemData,"R2");
            FieldNumbers[3] = mxGetFieldNumber(ElemData,"T1");
            FieldNumbers[4] = mxGetFieldNumber(ElemData,"T2");
			FieldNumbers[5] = mxGetFieldNumber(ElemData,"RApertures");
            FieldNumbers[6] = mxGetFieldNumber(ElemData,"EApertures");
            /* Fall through next section... */
            
	    case USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
                                 The second argument ponter to the array of field
                                 numbers is previously created with
                                 QuadLinPass( ..., MAKE_LOCAL_COPY)
                                 */
            le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
            
            /* Optional fields */
            pr1 = (FieldNumbers[1] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[1])) : NULL;
            pr2 = (FieldNumbers[2] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[2])) : NULL;
            pt1 = (FieldNumbers[3] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[3])) : NULL;
            pt2 = (FieldNumbers[4] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[4])) : NULL;
            RApertures = (FieldNumbers[5] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[5])) : NULL;
			EApertures = (FieldNumbers[6] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[6])) : NULL;
            break;
	    default:
            mexErrMsgTxt("No match for calling mode in function QuadMPoleFringePass\n");
	}
	DriftPass(r_in, le, pt1, pt2, pr1, pr2, RApertures, EApertures, num_particles);
	return FieldNumbers;
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        double *pr1, *pr2, *pt1, *pt2;
        double *RApertures, *EApertures;
        mxArray *tmpmxptr;

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
        DriftPass(r_in, le, pt1, pt2, pr1, pr2, RApertures, EApertures, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(1,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
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
