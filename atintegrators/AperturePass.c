
#include "atelem.c"
#include "atlalib.c"

struct elem {
    double *limits;
};

void AperturePass(double *r_in, double *limits, int num_particles)
{
    /* Checks X,Y of each input 6-vector and marks the corresponding element in
     * lossflag array with 0 if X,Y are exceed the limits given by limits
     * arrays
     * limits has 4 elements: (MinX, MaxX, MinY, MaxY)
     */
    int c;
    double *r6;
    for (c=0; c<num_particles; c++) {
        r6 = r_in+c*6;
        if (!atIsNaN(r6[0]))
        {
            /*  check if this particle is already marked as lost */
            checkiflostRectangularAp(r6, limits);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double *limits=atGetDoubleArray(ElemData,"Limits"); check_error();
	Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->limits=limits;
	}
    AperturePass(r_in, Elem->limits, num_particles);

    return Elem;
}

MODULE_DEF(AperturePass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double *limits=atGetDoubleArray(ElemData,"Limits"); check_error();
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        AperturePass(r_in,limits, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(1,1);
        mxSetCell(plhs[0],0,mxCreateString("Limits"));
        if(nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
