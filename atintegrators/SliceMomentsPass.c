#include "atelem.c"
#include "atimplib.c"
#include <math.h>
#include <float.h>
#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif


struct elem
{
  int startturn;
  int endturn;
  int nslice;
  double *stds;
  double *means;
  double *weight;
  double *z_cuts;
};


static void slice_bunch(double *r_in,int num_particles,int nslice,int turn,
                        int nturns, int nbunch, double *weight, double *means, double *stds,
                        double *z_cuts){

    int i,ii,ib;
    double *rtmp;

    double *smin = atMalloc(nbunch*sizeof(double));
    double *smax = atMalloc(nbunch*sizeof(double));
    double *hz = atMalloc(nbunch*sizeof(double));
    double *np_bunch = atMalloc(nbunch*sizeof(double));
    getbounds(r_in,nbunch,num_particles,smin,smax,z_cuts);

    for(i=0;i<nbunch;i++){
        hz[i] = (smax[i]-smin[i])/(nslice);
        np_bunch[i] = 0.0;
    }

    double *xpos = turnhistory + turn*nslice*nbunch;
    double *pxpos = turnhistory + (nturns + turn)*nslice*nbunch;
    double *ypos = turnhistory + (2*nturns + turn)*nslice*nbunch;
    double *pypos = turnhistory + (3*nturns + turn)*nslice*nbunch;
    double *dppos = turnhistory + (4*nturns + turn)*nslice*nbunch;
    double *zpos = turnhistory + (5*nturns + turn)*nslice*nbuncih;


    double *weight = turnhistory + (4*nturns-1)*nslice*nbunch;


    /*slices sorted from head to tail (increasing ct)*/
    for (i=0;i<num_particles;i++) {
        rtmp = r_in+i*6;
        ib = i%nbunch;
        np_bunch[ib] += 1.0;
        if (!atIsNaN(rtmp[0])) {
            register double x = rtmp[0];
            register double y = rtmp[2];
            register double ct = rtmp[5];
            if (ct < smin[ib]) {
                pslice[i] = ib*nslice;
            }
            else if (ct >smax[ib]){
                pslice[i] = nslice-1 + ib*nslice;
            }
            else if (ct == smax[ib]){
                ii = nslice-1 + ib*nslice;
                weight[ii] += 1.0;
                xpos[ii] += x;
                ypos[ii] += y;
                zpos[ii] += ct;
                pslice[i] = ii;
            }
            else {
                ii = (int)(floor((ct-smin[ib])/hz[ib])) + ib*nslice;
                weight[ii] += 1.0;
                xpos[ii] += x;
                ypos[ii] += y;
                zpos[ii] += ct;
                pslice[i] = ii;
            }
        }
    }

    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE,np_bunch,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,xpos,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ypos,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,zpos,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,weight,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    /*Compute average x/y position and weight of each slice */
    for (i=0;i<nslice*nbunch;i++) {
        ib = (int)(i/nslice);
        zpos[i] =  (weight[i]>0.0) ? zpos[i]/weight[i] : smin[ib]+(i%nslice+0.5)*hz[ib];
        zpos[i] += bunch_spos[ib]-bunch_spos[nbunch-1];
        xpos[i] =  (weight[i]>0.0) ? xpos[i]/weight[i] : 0.0;
        ypos[i] =  (weight[i]>0.0) ? ypos[i]/weight[i] : 0.0;
        weight[i] *= bunch_currents[ib]/np_bunch[ib];
    }
    atFree(np_bunch);
    atFree(smin);
    atFree(smax);
    atFree(hz);
};


void SliceMomentsPass(double *r_in, int nbunch, int num_particles, struct elem *Elem) {

    int startturn = Elem->startturn;
    int endturn = Elem->endturn;
    int nslice = Elem->nslice;
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