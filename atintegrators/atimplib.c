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
        ytmp0 = turnhistory + i*nslice+nslice*nturns;
        ytmp = turnhistory + (i+1)*nslice+nslice*nturns;
        ztmp0 = turnhistory + i*nslice+nslice*nturns*2;
        ztmp = turnhistory + (i+1)*nslice+nslice*nturns*2;
        wtmp0 = turnhistory + i*nslice+nslice*nturns*3;
        wtmp = turnhistory + (i+1)*nslice+nslice*nturns*3;
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
    ytmp = turnhistory + (nturns-1)*nslice+nslice*nturns;
    ztmp = turnhistory + (nturns-1)*nslice+nslice*nturns*2;
    wtmp = turnhistory + (nturns-1)*nslice+nslice*nturns*3; 
    for(ii=0;ii<nslice;ii++){
        xtmp[ii]=0.0;
        ytmp[ii]=0.0;
        ztmp[ii]=0.0;
        wtmp[ii]=0.0;
    }        
};



double *getbounds(double *r_in,int num_particles){
    double *rtmp;
    int i;
    static double bounds[2];
    double smin= DBL_MAX;
    double smax= -DBL_MAX;
    /*First find the min and the max of the distribution*/  
    for (i=0;i<num_particles;i++) {
        rtmp = r_in+i*6;
        if (!atIsNaN(rtmp[0])) {
            register double ct = rtmp[5];
            if (ct>smax) smax = ct;
            if (ct<smin) smin = ct;
        }
    }

    #ifdef MPI
    int flag;
    MPI_Initialized(&flag);
    if(flag){
        MPI_Allreduce(MPI_IN_PLACE,&smin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&smax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
        MPI_Barrier(MPI_COMM_WORLD);
    };
    #endif

    bounds[0]=smin;
    bounds[1]=smax;
    return bounds;
}


void slice_bunch(double *r_in,int num_particles,int nslice,int nturns,
                 double *turnhistory,int *pslice,double *z_cuts){
    
    int i,ii;
    double *rtmp;

    double *bounds = z_cuts ? z_cuts : getbounds(r_in,num_particles); 
    double smin = bounds[0];
    double smax = bounds[1];

    double hz = (smax-smin)/(nslice);
    double *xpos = turnhistory + (nturns-1)*nslice;
    double *ypos = turnhistory + (nturns-1)*nslice+nslice*nturns;
    double *zpos = turnhistory + (nturns-1)*nslice+nslice*nturns*2;
    double *weight = turnhistory + (nturns-1)*nslice+nslice*nturns*3;


    /*slices sorted from head to tail (increasing ct)*/
    for (i=0;i<num_particles;i++) {
        rtmp = r_in+i*6;
        if (!atIsNaN(rtmp[0])) {
            register double x = rtmp[0];
            register double y = rtmp[2];
            register double ct = rtmp[5];
            if (ct < smin) {
                pslice[i] = 0;
            }
            else if (ct >= smax) {
                pslice[i] = nslice-1;
            }
            else {
                ii = (int)floor((ct-smin)/hz);
                weight[ii] += 1.0;
                xpos[ii] += x;
                ypos[ii] += y;
                zpos[ii] += ct;
                pslice[i] = ii;                
            }
        }
    }


    double np_total = (double)num_particles;

    #ifdef MPI
    int flag;
    MPI_Initialized(&flag); 
    if(flag){
        int mpsize;
        MPI_Comm_size(MPI_COMM_WORLD,&mpsize); 
        np_total *= (double)mpsize;      
        MPI_Allreduce(MPI_IN_PLACE,xpos,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,ypos,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,zpos,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,weight,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    #endif

    /*Compute average x/y position and weight of each slice */
    for (i=0;i<nslice;i++) {
        zpos[i] =  (weight[i]>0.0) ? zpos[i]/weight[i] : (smin + (i+0.5)*hz);
        xpos[i] =  (weight[i]>0.0) ? xpos[i]/weight[i] : 0.0;
        ypos[i] =  (weight[i]>0.0) ? ypos[i]/weight[i] : 0.0;
        weight[i] /= np_total;
    } 

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
    #ifdef MPI
    int flag;
    MPI_Initialized(&flag); 
    if(flag){
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    }
    #endif
    for(i=nslice*(nturns-1);i<nslice*nturns;i++){  
        if(turnhistoryW[i]>0.0 && rank==(i+size)%size){
            for (ii=0;ii<nslice*nturns;ii++){
                ds = turnhistoryZ[ii]-turnhistoryZ[i];
                if(turnhistoryW[ii]>0.0 && -ds>=waketableT[0] && -ds<waketableT[nelem-1]){
                    wi = turnhistoryW[ii];
                    dx = turnhistoryX[ii];
                    dy = turnhistoryY[ii];
                    index = binarySearch(waketableT,-ds,nelem,0,0);              
                    if(waketableDX)kx[i-nslice*(nturns-1)] += dx*normfact[0]*wi*getTableWake(waketableDX,waketableT,-ds,index);
                    if(waketableDY)ky[i-nslice*(nturns-1)] += dy*normfact[1]*wi*getTableWake(waketableDY,waketableT,-ds,index);
                    if(waketableQX)kx2[i-nslice*(nturns-1)] += normfact[0]*wi*getTableWake(waketableQX,waketableT,-ds,index);
                    if(waketableQY)ky2[i-nslice*(nturns-1)] += normfact[1]*wi*getTableWake(waketableQY,waketableT,-ds,index);
                    if(waketableZ) kz[i-nslice*(nturns-1)] += normfact[2]*wi*getTableWake(waketableZ,waketableT,-ds,index);
                }            
            }
        }
    }
    #ifdef MPI
    if(flag){
        MPI_Allreduce(MPI_IN_PLACE,kx,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,ky,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,kx2,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,ky2,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,kz,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    #endif
};



