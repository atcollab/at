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


void rotate_table_history(long nturns,long nslice,double *turnhistoryX,double *turnhistoryY,double *turnhistoryZ,
                        double *turnhistoryW,double circumference){
    double *xtmp,*xtmp0;
    double *ytmp,*ytmp0;
    double *ztmp,*ztmp0;
    double *wtmp,*wtmp0;
    int i, ii;
    /*First rotate array*/
    for (i=0;i<nturns-1;i++){
        xtmp = turnhistoryX + (i+1)*nslice;
        xtmp0 = turnhistoryX + i*nslice;
        ytmp = turnhistoryY + (i+1)*nslice;
        ytmp0 = turnhistoryY + i*nslice;
        ztmp = turnhistoryZ + (i+1)*nslice;
        ztmp0 = turnhistoryZ + i*nslice;
        wtmp = turnhistoryW + (i+1)*nslice;
        wtmp0 = turnhistoryW + i*nslice;
        for(ii=0;ii<nslice;ii++){
            xtmp0[ii]=xtmp[ii];
            ytmp0[ii]=ytmp[ii];
            /*shift zpos by one turn*/
            ztmp0[ii]=ztmp[ii]+circumference;
            wtmp0[ii]=ytmp[ii];
        }
    }
    /*Now set last row to 0 (present slices)*/
    xtmp = turnhistoryX + (nturns-1)*nslice;
    ytmp = turnhistoryY + (nturns-1)*nslice;
    ztmp = turnhistoryZ + (nturns-1)*nslice;
    wtmp = turnhistoryW + (nturns-1)*nslice;
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
                 double *turnhistoryX,double *turnhistoryY,double *turnhistoryZ,
                 double *turnhistoryW,int *pslice,double *z_cuts){

    
    int i,ii;
    double *rtmp;
    double smin, smax;
    if (z_cuts) {
      smin=z_cuts[0];
      smax=z_cuts[1];
                 }
    else {
      double *bounds = getbounds(r_in,num_particles);
      smin = bounds[0];
      smax = bounds[1]; 
         }

    double hz;
    hz = (smax-smin)/(nslice);
    double *xpos = turnhistoryX + (nturns-1)*nslice;
    double *ypos = turnhistoryY + (nturns-1)*nslice;
    double *zpos = turnhistoryZ + (nturns-1)*nslice;
    double *weight = turnhistoryW + (nturns-1)*nslice;


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



