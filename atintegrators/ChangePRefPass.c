#include "atelem.c"
#include "atlalib.c"

struct elem {
    double Scaling;
};

void ChangePRefPass(double *r_in, double scaling, int num_particles)
/* scaling: scaling of the reference momentum
   r_in:    6-by-N matrix of initial conditions reshaped into 1-d array of 6*N elements
*/
{
    int c;

    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) \
    default(none) shared(r_in, num_particles, scaling ) private(c)
    for (c = 0; c<num_particles; c++) { /*Loop over particles  */
        double *r6 = r_in + 6*c;
        if (!atIsNaN(r6[0])) {
            ATChangePRef(r6, scaling);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Scaling;
        Scaling=atGetDouble(ElemData,"FieldScaling"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Scaling=Scaling;
    }
    ChangePRefPass(r_in, Elem->Scaling, num_particles);
    return Elem;
}

MODULE_DEF(ChangePRefPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Scaling;
        Scaling=atGetDouble(ElemData,"FieldScaling"); check_error();
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        ChangePRefPass(r_in, Scaling, num_particles);
    } else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(1,1);
        mxSetCell(plhs[0],0,mxCreateString("FieldScaling"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(0,1);
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
