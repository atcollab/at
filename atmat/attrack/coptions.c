#include <mex.h>
#ifndef OCTAVE
#include <matrix.h>
#endif

/*!
 * Return a structure describing the C compile options
 
 @param[out] options: structure with fields
                - "openmp": true if compiled for OpenMP

*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const char *fieldnames[] = {"openmp"};
    int omp;

    plhs[0] = mxCreateStructMatrix(1, 1, 1, fieldnames);
#ifdef _OPENMP
    omp = 1;
#else
    omp = 0;
#endif /*_OPENMP*/
    mxSetField(plhs[0], 0, "openmp", mxCreateDoubleScalar((double)omp));
}