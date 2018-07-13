/* IdentityPass.c 
   Accelerator Toolbox
   Revision 7/16/03
   A.Terebilo terebilo@slac.stanford.edu
*/

#include "atelem.c"
#include "atlalib.c"

struct elem 
{
  double *R1;
  double *R2;
  double *T1;
  double *T2;
  double *EApertures;
  double *RApertures;
};

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
        if (!atIsNaN(r6[0])) {
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

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
    }
    IdentityPass(r_in,Elem->T1,Elem->T2,Elem->R1,Elem->R2,
            Elem->RApertures,Elem->EApertures,num_particles);
    return Elem;
}

MODULE_DEF(IdentityPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        IdentityPass(r_in,T1,T2,R1,R2,RApertures,EApertures,num_particles);
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
#endif /*defined(MATLAB_MEX_FILE)*/
