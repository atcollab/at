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
  double *sizes;
  double *positions;
};


void BeamMomentsPass(double *r_in, int nbunch, int num_particles, struct elem *Elem) {

    int turn = Elem->turn;
    double *sizes = Elem->sizes;
    double *positions = Elem->positions;
        
    int i, ii, ib; 
    size_t sz = 2*nbunch*6*sizeof(double) + nbunch*sizeof(double);
    void *buffer = atMalloc(sz);
    double *dptr = (double *) buffer;
    
    double *avep=dptr; dptr += nbunch*6;
    double *sizep=dptr; dptr += nbunch*6;
    double *nparts=dptr;
    
    for (i=0; i<nbunch; i++) {
        nparts[i] = 0.0;
        for(ii=0; ii<6; ii++) {
            avep[i*6+ii] = 0.0;
            sizep[i*6+ii] = 0.0;
        }
    }
    
    for (i=0; i<num_particles; i++) {
        double *r6 = r_in+i*6;
        if(!atIsNaN(r6[0])){
            ib = i%nbunch;
            nparts[ib] += 1;
            for(ii=0; ii<6; ii++) {
                avep[6*ib+ii] += r6[ii];
                sizep[6*ib+ii] += r6[ii]*r6[ii];
            }
        }
    }   
    
    #ifdef MPI     
    MPI_Allreduce(MPI_IN_PLACE,sizep,6*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,avep,6*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,nparts,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
    for (i=0; i<nbunch; i++){       
        for (ii=0; ii<6; ii++){
            avep[6*i+ii] = avep[6*i+ii]/nparts[i];
            sizep[6*i+ii] = sqrt(sizep[6*i+ii]/nparts[i]-avep[6*i+ii]*avep[6*i+ii]); 
        }
    }    
     
    positions += 6*nbunch*turn;
    sizes += 6*nbunch*turn;
    memcpy(positions, avep, 6*nbunch*sizeof(double)); 
    memcpy(sizes, sizep, 6*nbunch*sizeof(double));     
}


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                                      double *r_in, int num_particles, struct parameters *Param)
{
    double *sizes;
    double *positions;  
    if (!Elem) {   
        positions=atGetDoubleArray(ElemData,"_positions"); check_error();
        sizes=atGetDoubleArray(ElemData,"_sizes"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->sizes=sizes;
        Elem->positions=positions;
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
    if (nrhs == 2) {
        double *sizes;
        double *positions; 
        positions=atGetDoubleArray(ElemData,"_positions"); check_error();
        sizes=atGetDoubleArray(ElemData,"_sizes"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->sizes=sizes;
        Elem->positions=positions;
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
        mxSetCell(plhs[0],0,mxCreateString("_positions"));
        mxSetCell(plhs[0],1,mxCreateString("_sizes"));
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif
