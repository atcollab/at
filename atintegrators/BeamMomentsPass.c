#include "atelem.c"
#include "atimplib.c"
#include <math.h>
#include <float.h>
#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif
/*
 * Beam moments pass method by Simon White.  
 *
 */

struct elem
{
  int turn;
  double *stds;
  double *means;
};


void BeamMomentsPass(double *r_in, int nbunch, int num_particles, struct elem *Elem) {

    int turn = Elem->turn;
    double *stds = Elem->stds;
    double *means = Elem->means;
        
    int i, ii, ib; 
    void *buffer = atCalloc(2*nbunch*6+nbunch, sizeof(double));
    double *dptr = (double *) buffer;
    
    double *meanp=dptr; dptr += nbunch*6;
    double *stdp=dptr; dptr += nbunch*6;
    double *nparts=dptr;
    
    for (i=0; i<num_particles; i++) {
        double *r6 = r_in+i*6;
        if(!atIsNaN(r6[0])){
            ib = i%nbunch;
            nparts[ib] += 1;
            for(ii=0; ii<6; ii++) {
                meanp[6*ib+ii] += r6[ii];
                stdp[6*ib+ii] += r6[ii]*r6[ii];
            }
        }
    }   
    
    #ifdef MPI     
    MPI_Allreduce(MPI_IN_PLACE,stdp,6*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,meanp,6*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,nparts,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    for (i=0; i<nbunch; i++){
        for (ii=0; ii<6; ii++){
            meanp[6*i+ii] = meanp[6*i+ii]/nparts[i];
            stdp[6*i+ii] = sqrt(stdp[6*i+ii]/nparts[i]-meanp[6*i+ii]*meanp[6*i+ii]); 
        }
    }    

    means += 6*nbunch*turn;
    stds += 6*nbunch*turn;
    memcpy(means, meanp, 6*nbunch*sizeof(double));
    memcpy(stds, stdp, 6*nbunch*sizeof(double));
    atFree(buffer);  
}


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                                      double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double *means;
        double *stds;
        int ndim = 3;
        int dims[] = {6, Param->nbunch, Param->num_turns};
        means=atGetDoubleArray(ElemData,"_means"); check_error();
        stds=atGetDoubleArray(ElemData,"_stds"); check_error();
        atCheckArrayDims(ElemData,"_means", ndim, dims); check_error();
        atCheckArrayDims(ElemData,"_stds", ndim, dims); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->stds=stds;
        Elem->means=means;
        Elem->turn = 0;
    }
    BeamMomentsPass(r_in, Param->nbunch, num_particles, Elem);
    Elem->turn++;
    return Elem;
}

MODULE_DEF(BeamMomentsPass)        /* Dummy module initialisation */
#endif

#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
    
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        struct elem El, *Elem=&El;
        
        double *means;
        double *stds; 
        means=atGetDoubleArray(ElemData,"_means"); check_error();
        stds=atGetDoubleArray(ElemData,"_stds"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->stds=stds;
        Elem->means=means;
        Elem->turn = 0;
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix: particle array");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        BeamMomentsPass(r_in,1,num_particles,Elem);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(2,1);
        mxSetCell(plhs[0],0,mxCreateString("_means"));
        mxSetCell(plhs[0],1,mxCreateString("_stds"));
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif
