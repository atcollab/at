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


void slice_bunch(double *r_in,int num_particles,int nslice,double *bounds,double *weight,
                 double *xpos,double *ypos,double *zpos,int *countslc,int *pslice){

    
    int i,ii;
    double *rtmp;
    double smin = bounds[0];
    double smax = bounds[1];  
    double hz;
    hz = (smax-smin)/(nslice);


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
                countslc[ii]++;
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
        MPI_Allreduce(MPI_IN_PLACE,countslc,nslice,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    #endif

    /*Compute average x/y position and weight of each slice */
    for (i=0;i<nslice;i++) {
        double count = (double) countslc[i];
        weight[i]=count/np_total;
        zpos[i] =  (countslc[i]>0) ? zpos[i]/count : (smin + (i+0.5)*hz);
        xpos[i] =  (countslc[i]>0) ? xpos[i]/count : 0.0;
        ypos[i] =  (countslc[i]>0) ? ypos[i]/count : 0.0;
    } 

};



