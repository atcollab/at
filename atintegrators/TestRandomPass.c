 /* RandomPass.c
   Accelerator Toolbox

   Test of random generators
 */

#include "atelem.c"
#include "atlalib.c"
#include "atrandom.c"
#ifdef MPI
#include <mpi.h>
#endif

struct elem 
{
    int dummy;  /* to make Windows compiler happy */
};

static void RandomPass(double *r_in,
        pcg32_random_t* common_rng,
        pcg32_random_t* thread_rng,
        int num_particles)
{	
    double common_val;
#ifdef MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    int rank = 0;
#endif /* MPI */

    common_val = atrandn_r(common_rng, 0.0, 0.001);
    for (int c = 0; c<num_particles; c++) {	/*Loop over particles  */
        double *r6 = r_in+c*6;
        r6[0] = atrandn_r(thread_rng, 0.0, 0.001);
        r6[2] = common_val;
        r6[4] = 0.01*rank;
        r6[5] = 0.01*c;
    }

    common_val = atrandn_r(common_rng, 0.0, 0.001);
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none)                      \
    shared(r_in, num_particles, common_val, thread_rng)
    for (int c = 0; c<num_particles; c++) {	/*Loop over particles  */
        double *r6 = r_in+c*6;
        r6[1] = atrandn_r(thread_rng, 0.0, 0.001);;
        r6[3] = common_val;
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
    }
    RandomPass(r_in,Param->common_rng,Param->thread_rng,num_particles);
    return Elem;
}

MODULE_DEF(TestRandomPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        RandomPass(r_in,&pcg32_global,&pcg32_global,num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(1,0);
        
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(1,0);
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
