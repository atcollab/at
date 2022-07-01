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
 * Beam monitor pass method by Simon White.  
 * User may contact simon.white@esrf.fr for questions and comments.
 */

struct elem
{
  int nbunch;
  int turn;
  double *sizes;
  double *positions;
};


void BeamMonitorPass(double *r_in, int num_particles, struct elem *Elem) {


    long nbunch = Elem->nbunch;
    int turn = Elem->turn;
    double *sizes = Elem->sizes;
    double *positions = Elem->positions;
        
    int i, ii, ib; 
    double nparts[nbunch];
    double avep[nbunch*6];
    double sizep[nbunch*6];
    
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
    long nbunch, nturns;
    double *sizes;
    double *positions;  
    if (!Elem) {   
        positions=atGetDoubleArray(ElemData,"_positions"); check_error();
        sizes=atGetDoubleArray(ElemData,"_sizes"); check_error();
        nbunch=atGetLong(ElemData,"_nbunch"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->nbunch=nbunch;
        Elem->sizes=sizes;
        Elem->positions=positions;
        Elem->turn = 0;
    }
    BeamMonitorPass(r_in, num_particles, Elem);
    Elem->turn++;
    return Elem;
}

MODULE_DEF(BeamMonitorPass)        /* Dummy module initialisation */
#endif
