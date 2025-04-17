#include "atconstants.h"
#include "atelem.c"
#include <math.h>
#include <float.h>
#include <complex.h>
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


static double getTableWake(double *waketable,double *waketableT,double distance,int index){
    double w = waketable[index] + (distance-waketableT[index])*(waketable[index+1]-waketable[index])/
          (waketableT[index+1]-waketableT[index]);
    if(atIsNaN(w)){
        return 0;
    }else{
        return w;
    };
};

static void rotate_table_history(long nturns,long nslice,double *turnhistory,double circumference){

    double *xtmp,*xtmp0;
    double *ytmp,*ytmp0;
    double *ztmp,*ztmp0;
    double *wtmp,*wtmp0;
    double *t0;
    int i, ii;    
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
            ztmp0[ii]=ztmp[ii]-circumference;
            wtmp0[ii]=wtmp[ii];
        }
    }
    
    for(ii=1;ii<5; ii++){ 
        t0 = turnhistory + (ii*nturns-1)*nslice;
        for(i=0; i<nslice; i++){
            t0[i] = 0.0;
        }
    }
};

static void getbounds(double *r_in, int nbunch, int num_particles, double *smin,
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


static void slice_bunch(double *r_in,int num_particles,int nslice,int nturns,
                 int nbunch,double *bunch_spos,double *bunch_currents,
                 double *turnhistory,int *pslice,double *z_cuts){
    
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
        if (np_bunch[ib] == 0.0) {
            weight[i] = 0.0;
            }
        else {
            weight[i] *= bunch_currents[ib]/np_bunch[ib];
        }
    } 
    atFree(np_bunch);
    atFree(smin);
    atFree(smax);
    atFree(hz);
};

static void compute_kicks(int nslice,int nturns,int nelem,
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
    if(waketableDX)MPI_Allreduce(MPI_IN_PLACE,kx,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(waketableDY)MPI_Allreduce(MPI_IN_PLACE,ky,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(waketableQX)MPI_Allreduce(MPI_IN_PLACE,kx2,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(waketableQY)MPI_Allreduce(MPI_IN_PLACE,ky2,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(waketableZ)MPI_Allreduce(MPI_IN_PLACE,kz,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
};


static void wakefunc_long_resonator(double ds, double freqres, double qfactor, double rshunt, double beta, double *wake) {

    double omega, alpha, omegabar;
    wake[0] = 0.0;
    wake[1] = 0.0;
    double dt;
    
    omega = TWOPI * freqres;
    alpha = omega / (2 * qfactor);
    omegabar = sqrt(fabs(omega*omega - alpha*alpha));
    
    dt = -ds/(beta * C0);      
             
    if (dt==0) {
        wake[0] = rshunt * alpha;
        wake[1] = 0.0;
    } else if (dt<0) {
        if (qfactor > 0.5) {
            wake[0] = 2 * rshunt * alpha * exp(alpha * dt) * (cos(omegabar * dt) + \
                      alpha / omegabar * sin(omegabar*dt));
            wake[1] = 2 * rshunt * alpha * exp(alpha * dt) * (sin(omegabar * dt) - \
                      alpha / omegabar * cos(omegabar*dt));
        } else if (qfactor == 0.5) {
            wake[0] = 2 * rshunt * alpha * exp(alpha * dt) * (1. + alpha * dt);
            wake[1] = 0.0;
        } else if (qfactor < 0.5) {
            wake[0] = 2 * rshunt * alpha * exp(alpha * dt) * (cosh(omegabar * dt) + \
                      alpha / omegabar * sinh(omegabar * dt)); 
            wake[1] = 2 * rshunt * alpha * exp(alpha * dt) * (sinh(omegabar * dt) - \
                      alpha / omegabar * cosh(omegabar*dt));
        }       
    } else {
        wake[0] = 0.0;
        wake[1] = 0.0;
    }            
}


static void compute_kicks_longres(int nslice,int nbunch,int nturns, double *turnhistory,double normfact,
                           double *kz,double freq, double qfactor, double rshunt,
                           double beta, double *vbeamk, double energy, double *vbunch) {

    int rank=0;
    int size=1;
    int i,ii,ib;
    double ds,wi,wii;
    double *turnhistoryZ = turnhistory+nslice*nbunch*nturns*2;
    double *turnhistoryW = turnhistory+nslice*nbunch*nturns*3;
    double wake[2];
    double vba, vbp;
    double *vbr = vbunch;
    double *vbi = vbunch+nbunch;
    double totalW = 0;
    double *totalWb = atMalloc(nbunch*sizeof(double));
    
    for (i=0;i<nslice*nbunch;i++) {
        ib = (int)(i/nslice);
        kz[i]=0.0;
        vbr[ib] = 0.0;
        vbi[ib] = 0.0;
        totalWb[ib] = 0.0;
    }

    vbeamk[0] = 0.0;
    vbeamk[1] = 0.0;
    wake[0] = 0.0;
    wake[1] = 0.0;


    #ifdef MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    #endif
    for(i=nslice*nbunch*(nturns-1);i<nslice*nbunch*nturns;i++){  
        ib = (int)((i-nslice*nbunch*(nturns-1))/nslice);
        wi = turnhistoryW[i];
        if(turnhistoryW[i]>0.0 && rank==(i+size)%size){
            totalW += wi;
            totalWb[ib] += wi;
            for (ii=0;ii<nslice*nbunch*nturns;ii++){
                ds = turnhistoryZ[i]-turnhistoryZ[ii];
                if(turnhistoryW[ii]>0.0 && ds>=0){
                    wii = turnhistoryW[ii];
                    wakefunc_long_resonator(ds,freq,qfactor,rshunt,beta,wake);       
                    kz[i-nslice*nbunch*(nturns-1)] += normfact*wii*wake[0];
                    vbeamk[0] += normfact*wii*wake[0]*energy*wi;
                    vbeamk[1] -= normfact*wii*wake[1]*energy*wi;
                    vbr[ib] += normfact*wii*wake[0]*energy*wi;
                    vbi[ib] -= normfact*wii*wake[1]*energy*wi;
                }            
            }
        }
    }

    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE,kz,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,vbeamk,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,vbr,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,vbi,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&totalW,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,totalWb,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
    vba = sqrt(vbeamk[0]*vbeamk[0]+vbeamk[1]*vbeamk[1])/totalW;
    vbp = atan2(vbeamk[1],vbeamk[0]);
    vbeamk[0] = vba;
    vbeamk[1] = vbp;
    
    for(i=0;i<nbunch;i++){
        double vr = vbr[i]/totalWb[i];
        double vi = vbi[i]/totalWb[i];
        vbr[i] = sqrt(vr*vr+vi*vi); 
        vbi[i] = atan2(vi,vr);
    }
    atFree(totalWb);
};


static void compute_kicks_phasor(int nslice, int nbunch, int nturns, double *turnhistory,
                          double normfact, double *kz,double freq, double qfactor,
                          double rshunt, double *vbeam, double circumference,
                          double energy, double beta, double *vbeamk, double *vbunch){ 
                          
    #ifndef _MSC_VER  
    int i,ib,is;
    double wi;
    double selfkick;
    int sliceperturn = nslice*nbunch;
    double dt =0.0;
    double *turnhistoryZ = turnhistory+nslice*nbunch*nturns*2;
    double *turnhistoryW = turnhistory+nslice*nbunch*nturns*3;
    double omr = TWOPI*freq;
    double complex vbeamc = vbeam[0]*cexp(_Complex_I*vbeam[1]);
    double complex vbeamkc = 0.0;
    double kloss = rshunt*omr/(2*qfactor);
    double bc = beta*C0;
    double *vbr = vbunch;
    double *vbi = vbunch+nbunch;
    double totalW=0.0;
    double *totalWb = atMalloc(nbunch*sizeof(double));
    
    for (i=0;i<sliceperturn;i++) {
        ib = (int)(i/nslice);
        kz[i]=0.0;
        vbr[ib] = 0.0;
        vbi[ib] = 0.0;
        totalWb[ib] = 0.0;
    }
    
    for(i=sliceperturn*(nturns-1);i<sliceperturn*nturns;i++){
        ib = (int)((i-sliceperturn*(nturns-1))/nslice);
        wi = turnhistoryW[i];
        selfkick = normfact*wi*kloss*energy;
        if(i==sliceperturn*(nturns-1)){
            /*At the end of the turn, the vbeamc is
            reverted to -final value, which stores the
            dt information from previous turn. This extra
            circumference is needed to take this into account. */
            dt = (circumference+turnhistoryZ[i])/bc;
        }else{
            /* This is dt between each slice*/
            dt = (turnhistoryZ[i]-turnhistoryZ[i-1])/bc;
        }
        vbeamc *= cexp((_Complex_I*omr-omr/(2*qfactor))*dt);
        /*vbeamkc is average kick i.e. average vbeam*/   
        vbeamkc += (vbeamc+selfkick)*wi;
        totalW += wi;
        totalWb[ib] += wi;
        kz[i-sliceperturn*(nturns-1)] = creal((vbeamc + selfkick)/energy);
        vbr[ib] += creal((vbeamc + selfkick)*wi);
        vbi[ib] += cimag((vbeamc + selfkick)*wi);
        vbeamc += 2*selfkick;    
    }
    
    /*This takes the vbeam backwards in time to effectively store the
    final slice position */
    dt = -turnhistoryZ[sliceperturn*nturns-1]/bc;    
    vbeamc *= cexp((_Complex_I*omr-omr/(2*qfactor))*dt);

    vbeam[0] = cabs(vbeamc);
    vbeam[1] = carg(vbeamc);
    vbeamkc /= (totalW);
    vbeamk[0] = cabs(vbeamkc);
    vbeamk[1] = carg(vbeamkc);   
    
    for(i=0;i<nbunch;i++){
        double vr = vbr[i]/totalWb[i];
        double vi = vbi[i]/totalWb[i];
        vbr[i] = sqrt(vr*vr+vi*vi); 
        vbi[i] = atan2(vi,vr);
    }
    atFree(totalWb);
    #endif    
};


static void update_vgen(double *vbeam,double *vcav,double *vgen,double voltgain,double phasegain,double detune_angle){
    double vbeamr = vbeam[0]*cos(vbeam[1]);
    double vbeami = vbeam[0]*sin(vbeam[1]);
    double vcavr = vcav[0]*cos(vcav[1]);
    double vcavi = vcav[0]*sin(vcav[1]);   
    double vgenr = vcavr - vbeamr;
    double vgeni = vcavi - vbeami; 
    double vga = sqrt(vgenr*vgenr+vgeni*vgeni);   
    double vgp = atan2(vgeni,vgenr)-vcav[1]+detune_angle;
    vgen[0] += (vga-vgen[0])*voltgain;
    vgen[1] += (vgp-vgen[1])*phasegain;
}
