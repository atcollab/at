#include "atelem.c"
#include "atrandom.c"

struct elem {
    double* Lmatp;
};

void QuantDiffPass(double* r_in, double* Lmatp, int nturn,
    pcg32_random_t* rng,
    int num_particles)
    /* Lmatp 6x6 matrix
     * r_in - 6-by-N matrix of initial conditions reshaped into
     * 1-d array of 6*N elements
     */
{
/*  The behaviour of random generators with OpenMP is doubtful. OpenMP disabled until
    it's understood
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
    shared(r_in, num_particles, Lmatp, rng)
*/
    for (int c = 0; c < num_particles; c++) {
        /*Loop over particles  */
        int i, j;
        double randnorm[6];
        double diffusion[6];
        double* r6 = r_in + c * 6;
        for (i = 0; i < 6; i++) {
            diffusion[i] = 0.0;
            randnorm[i] = atrandn_r(rng, 0.0, 1.0);
        }

        for (i = 0; i < 6; i++) {
            for (j = 0; j <= i; j++) {
                diffusion[i] += randnorm[j] * Lmatp[i + 6 * j];
            }
        }
        if (!atIsNaN(r6[0])) {
            r6[0] += diffusion[0];
            r6[1] += diffusion[1];
            r6[2] += diffusion[2];
            r6[3] += diffusion[3];
            r6[4] += diffusion[4];
            r6[5] += diffusion[5];
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem* trackFunction(const atElem* ElemData, struct elem* Elem,
    double* r_in, int num_particles, struct parameters* Param)
{
    int nturn = Param->nturn;
    if (!Elem) {
        double* Lmatp;
        Lmatp=atGetDoubleArray(ElemData,"Lmatp"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Lmatp = Lmatp;
    }
    QuantDiffPass(r_in, Elem->Lmatp, nturn, Param->thread_rng, num_particles);
    return Elem;
}

MODULE_DEF(QuantDiffPass) /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs >= 2) {
        double* r_in;
        const mxArray* ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double* Lmatp;
        Lmatp=atGetDoubleArray(ElemData,"Lmatp"); check_error();
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        QuantDiffPass(r_in, Lmatp, 0, &pcg32_global, num_particles);
    } else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(1, 1);
        mxSetCell(plhs[0], 0, mxCreateString("Lmatp"));
        if (nlhs > 1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(0, 1);
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
