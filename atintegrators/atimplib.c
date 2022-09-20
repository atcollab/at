#include "atelem.c"
#include <math.h>
#include <float.h>
#ifdef MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif


int binarySearch(double *array,double value,int upper,int lower,int nStep){
    int pivot = (int)(lower+upper)/2;
    if ((upper-lower)<=1){
        return lower;
    };
    if (value < array[pivot]){
        upper = pivot;
        nStep+=1;
        return binarySearch(array,value,upper,lower,nStep);
    } else if (value > array[pivot]){
        lower = pivot;
        nStep+=1;
        return binarySearch(array,value,upper,lower,nStep);
    }else{
        return pivot;
    };
};


double getTableWake(double *waketable,double *waketableT,double distance,int index){
    double w = waketable[index] + (distance-waketableT[index])*(waketable[index+1]-waketable[index])/
          (waketableT[index+1]-waketableT[index]);
    if(atIsNaN(w)){
        return 0;
    }else{
        return w;
    };
};


void rotate_table_history(long nturns,long nslice,double *turnhistory,double circumference){
    double *xtmp,*xtmp0;
    double *ytmp,*ytmp0;
    double *ztmp,*ztmp0;
    double *wtmp,*wtmp0;
    int i, ii;    
    /*First rotate array*/

    for (i=0;i<nturns-1;i++){
        xtmp0 = turnhistory + i*nslice;
        xtmp = turnhistory + (i+1)*nslice;
        ytmp0 = turnhistory + (i+nturns)*nslice;
        ytmp = turnhistory + (i+nturns+1)*nslice;
        ztmp0 = turnhistory + (i+2*nturns)*nslice;
        ztmp = turnhistory + (i+2*nturns+1)*nslice;
        wtmp0 = turnhistory + (i+3*nturns)*nslice;
        wtmp = turnhistory + (i+3*nturns+1)*nslice;
        for(ii=0;ii<nslice;ii++){
            xtmp0[ii]=xtmp[ii];
            ytmp0[ii]=ytmp[ii];
            /*shift zpos by one turn*/
            ztmp0[ii]=ztmp[ii]-circumference;
            wtmp0[ii]=wtmp[ii];
        }
    }
    /*Now set last row to 0 (present slices)*/
    
    xtmp = turnhistory + (nturns-1)*nslice;
    ytmp = turnhistory + (2*nturns-1)*nslice;
    ztmp = turnhistory + (3*nturns-1)*nslice;
    wtmp = turnhistory + (4*nturns-1)*nslice; 
    for(ii=0;ii<nslice;ii++){
        xtmp[ii]=0.0;
        ytmp[ii]=0.0;
        ztmp[ii]=0.0;
        wtmp[ii]=0.0;
    }        
};


void getbounds(double *r_in, int nbunch, int num_particles, double *smin,
               double *smax, double *z_cuts){
    double *rtmp;
    int i, ib;
    if(z_cuts){
        for(i=0;i<nbunch; i++){
            smin[i] = z_cuts[0];
            smax[i] = z_cuts[1];
        }
    }else{
        for(i=0;i<nbunch; i++){
            smin[i] = DBL_MAX;
            smax[i] = -DBL_MAX;
        }
        /*First find the min and the max of the distribution*/  
        for (i=0;i<num_particles;i++) {
            rtmp = r_in+i*6;
            ib = i%nbunch;
            if (!atIsNaN(rtmp[0])) {
                register double ct = rtmp[5];
                if (ct>smax[ib]) smax[ib] = ct;
                if (ct<smin[ib]) smin[ib] = ct;
            }
        }

        #ifdef MPI
        MPI_Allreduce(MPI_IN_PLACE,smin,nbunch,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,smax,nbunch,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
        MPI_Barrier(MPI_COMM_WORLD);
        #endif

        for(i=0;i<nbunch;i++){
            if(smin[i]==smax[i]){
                smin[i] -= 1.0e-12;
                smax[i] += 1.0e-12;
            }
        }
    }
}


void slice_bunch(double *r_in,int num_particles,int nslice,int nturns,
                 int nbunch,double *bunch_spos,double *bunch_currents,
                 double *turnhistory,int *pslice,double *z_cuts){
    
    int i,ii,ib,is;
    double *rtmp;
    
    double *smin = malloc(nbunch*sizeof(double));
    double *smax = malloc(nbunch*sizeof(double));
    double *hz = malloc(nbunch*sizeof(double));
    double *np_bunch = malloc(nbunch*sizeof(double));
    getbounds(r_in,nbunch,num_particles,smin,smax,z_cuts);     
    
    for(i=0;i<nbunch;i++){
        hz[i] = (smax[i]-smin[i])/(nslice);
    }

    double *xpos = turnhistory + (nturns-1)*nslice*nbunch;
    double *ypos = turnhistory + (2*nturns-1)*nslice*nbunch;
    double *zpos = turnhistory + (3*nturns-1)*nslice*nbunch;
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
    MPI_Allreduce(MPI_IN_PLACE,xpos,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ypos,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,zpos,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,weight,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
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
    free(np_bunch);
    free(smin);
    free(smax);
    free(hz);
};

void compute_kicks(int nslice,int nturns,int nelem,
                   double *turnhistory,double *waketableT,double *waketableDX,
                   double *waketableDY,double *waketableQX,double *waketableQY,
                   double *waketableZ,double *normfact, double *kx,double *ky,
                   double *kx2,double *ky2,double *kz){
    int rank=0;
    int size=1;
    int i,ii,index;
    double ds,wi,dx,dy;
    double *turnhistoryX = turnhistory;
    double *turnhistoryY = turnhistory+nslice*nturns;
    double *turnhistoryZ = turnhistory+nslice*nturns*2;
    double *turnhistoryW = turnhistory+nslice*nturns*3;

    for (i=0;i<nslice;i++) {
        kx[i]=0.0;
        ky[i]=0.0;
        kx2[i]=0.0;
        ky2[i]=0.0;
        kz[i]=0.0;
    }

    #ifdef MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    #endif
    for(i=nslice*(nturns-1);i<nslice*nturns;i++){  
        if(turnhistoryW[i]>0.0 && rank==(i+size)%size){
            for (ii=0;ii<nslice*nturns;ii++){
                ds = turnhistoryZ[i]-turnhistoryZ[ii];
                wi = turnhistoryW[ii];
                if(wi>0.0 && ds>=waketableT[0] && ds<waketableT[nelem-1]){
                    dx = turnhistoryX[ii];
                    dy = turnhistoryY[ii];
                    index = binarySearch(waketableT,ds,nelem,0,0);          
                    if(waketableDX)kx[i-nslice*(nturns-1)] += dx*normfact[0]*wi*getTableWake(waketableDX,waketableT,ds,index);
                    if(waketableDY)ky[i-nslice*(nturns-1)] += dy*normfact[1]*wi*getTableWake(waketableDY,waketableT,ds,index);
                    if(waketableQX)kx2[i-nslice*(nturns-1)] += normfact[0]*wi*getTableWake(waketableQX,waketableT,ds,index);
                    if(waketableQY)ky2[i-nslice*(nturns-1)] += normfact[1]*wi*getTableWake(waketableQY,waketableT,ds,index);
                    if(waketableZ) kz[i-nslice*(nturns-1)] += normfact[2]*wi*getTableWake(waketableZ,waketableT,ds,index);
                }            
            }
        }
    }
    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE,kx,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ky,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,kx2,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ky2,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,kz,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
};



