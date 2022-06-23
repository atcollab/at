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
  int nbunch
  double *sizes
  double *positions  
};


void BeamMonitorPass(double *r_in, int num_particles, struct elem *Elem) {


    long nbunch = Elem->nbunch;
    double *sizes = Elem->sizes;
    double *positions = Elem->positions;
        
    int i, ii, ib; 
    double *nparts = [nbunch];

    for (i=0; i<nbunch; i++) {
        npart[i] = 0.0;
    }
    
    for (i=0; i<num_particles; i++) {
        double *r6 = r_in+i*6;
        if(!atIsNaN(r6[0])){
            ib = i%nbunch;
            nparts[ib] += 1;
            for(ii=0; ii<6; ii++) {
                positions[6*ib+ii] += r6[ii]
                sizes[6*ib+ii] += r6[ii]*r6[ii]
            }
        }
    }   
    
    #ifdef MPI     
    MPI_Allreduce(MPI_IN_PLACE,sizes,6*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,positions,6*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,nparts,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
    for (i=0; i<nbunch; i++){       
        for (ii=0; ii<6; ii++){
            positions[6*i+ii] = positions[6*i+ii]/nparts[i];
            sizes[6*i+ii] = sqrt(sizes[6*i+ii]/nparts[i]-positions[6*i+ii]*positions[6*i+ii]); 
        }         
    }
        
};
