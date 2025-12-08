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
                   double *waketableZ, double *waketableCX, double *waketableCY,
                   double *normfact, double *kx,double *ky,
                   double *kx2,double *ky2,double *kz, double *kcx, double *kcy){
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
        kcx[i]=0.0;
        kcy[i]=0.0;
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
                    if(waketableCX)kcx[i-nslice*(nturns-1)] += normfact[0]*wi*getTableWake(waketableCX,waketableT,ds,index);
                    if(waketableCY)kcy[i-nslice*(nturns-1)] += normfact[1]*wi*getTableWake(waketableCY,waketableT,ds,index);
                    
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
    if(waketableCX)MPI_Allreduce(MPI_IN_PLACE,kcx,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(waketableCY)MPI_Allreduce(MPI_IN_PLACE,kcy,nslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
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


static void compute_kicks_phasor(int nslice, int nbunch, int nturns, double *turnhistory,
                          double normfact, double *vbeam_kicks, double resfreq, double qfactor,
                          double rshunt, double *vbeam_phasor, double circumference,
                          double energy, double beta, double *ave_vbeam, double *vbunch,
                          double *bunch_spos, int M, 
                          double *fillpattern){ 
                          
    #ifndef _MSC_VER  
    int i,ib;
    double wi;
    double selfkick;
    double dt =0.0;
    double *turnhistoryZ = turnhistory+nslice*nbunch*nturns*2;
    double *turnhistoryW = turnhistory+nslice*nbunch*nturns*3;
    double omr = TWOPI*resfreq;
    double complex vbeam_complex = vbeam_phasor[0]*cexp(_Complex_I*vbeam_phasor[1]);
    double kloss = rshunt*omr/(2*qfactor);
    double bc = beta*C0;
    double *vbr = vbunch;
    double *vbi = vbunch+nbunch;
    int ibunch, islice, total_slice_counter;
    int bunch_counter = 0;
    double bucket_curr = 0.0;
    double main_bucket = circumference / M;

    
    
    for (i=0;i<nslice*nbunch;i++) {
        ibunch = (int)(i/nslice);
        vbeam_kicks[i] = 0.0;
        vbr[ibunch] = 0.0;
        vbi[ibunch] = 0.0;
    }

    /* The vbeam_complex will always be sent to the center of the next bucket */

    for(ibunch=0; ibunch<M; ibunch++){
        bucket_curr = fillpattern[ibunch];
        
        if(bucket_curr!=0.0){
            for(islice=0; islice<nslice; islice++){
                total_slice_counter = islice + nslice*bunch_counter; 
                wi = turnhistoryW[total_slice_counter];
                selfkick = normfact*wi*kloss*energy; /*normfact*energy is -t0 . This number comes out to be negative, which is correct*/       
                if(islice==0){
                    /* TurnhistoryZ goes from -bucket991 to bucket0 */
                    dt = (turnhistoryZ[total_slice_counter] + bunch_spos[nbunch-1-bunch_counter])/bc;
                }else{
                    /* This is dt between each slice*/
                    dt = (turnhistoryZ[total_slice_counter]-turnhistoryZ[total_slice_counter-1])/bc;
                }
                
                /* track the dt */
                vbeam_complex *= cexp((_Complex_I*omr-omr/(2*qfactor))*dt);
                vbeam_kicks[total_slice_counter] = creal((vbeam_complex + selfkick)/energy);
                
                vbeam_complex += 2*selfkick;    
               
            }
            /* back to the center of the bucket */
            dt = -(turnhistoryZ[total_slice_counter] + bunch_spos[nbunch - 1 - bunch_counter])/bc;
            vbeam_complex *= cexp((_Complex_I*omr-omr/(2*qfactor))*dt);

            vbr[bunch_counter] = cabs(vbeam_complex);
            vbi[bunch_counter] = carg(vbeam_complex);

            bunch_counter += 1;
        }

        ave_vbeam[0] += cabs(vbeam_complex)/M;
        ave_vbeam[1] += carg(vbeam_complex)/M;
        
        /* advance the phasor to the center of the next bucket */
       
        dt = main_bucket/bc;
        vbeam_complex *= cexp((_Complex_I*omr-omr/(2*qfactor))*dt);
        
    }
    /* store the phasor for the next turn */
    vbeam_phasor[0] = cabs(vbeam_complex);
    vbeam_phasor[1] = carg(vbeam_complex);
    

    #endif    
};


static void update_vgen(double *vbeam,double *vcav,double *vgen, double voltgain,double phasegain,double detune_angle){

    double vbeamr_meas = vbeam[0]*cos(vbeam[1]);
    double vbeami_meas = vbeam[0]*sin(vbeam[1]);
    
    double vgenr_meas = -vgen[0]*sin(vgen[1]);
    double vgeni_meas = vgen[0]*cos(vgen[1]);      
    
    double vcavr_meas = vgenr_meas + vbeamr_meas;
    double vcavi_meas = vgeni_meas + vbeami_meas;   

    double vcav_meas = sqrt(vcavr_meas*vcavr_meas + vcavi_meas*vcavi_meas); 
    double phis_meas = -atan2(vcavr_meas, vcavi_meas);
        
    double phis = vcav[1];   
    double ptmp = phis_meas - phis; /* this applies to thetag*/
    
    double dttmp = vgen[1] - vgen[2] - phis + detune_angle;

    double dtmp = vcav[0] / vcav_meas;
    /*printf("vcav_amp_set: %f, vcav_meas: %f \n", vcav[0], vcav_meas);*/
    
    vgen[3] *= pow(dtmp,voltgain);
    vgen[2] += dttmp*phasegain; 
    vgen[1] -= ptmp*phasegain;
    vgen[0] = vgen[3]*cos(vgen[2]);
}


static void update_passive_frequency(double *vbeam, double *vcav, double *vgen, double phasegain){
    /* The cavity voltage is
    V(t) = 2*I0*rs*cos(psi)*exp(i(wt+psi))
    We save the amplitude of vbeam, so the exponent goes to 1.
    Therefore vbeam[0] = 2*I0*rs*cos(psi) which is the cavity voltage.
    */
    double vset = vcav[0];
    double psi = vgen[2];
    double vpeak = vbeam[0]; /* Peak amplitude of cavity voltage */
    double delta_v = vset - vpeak;
    double grad = vbeam[0]*sin(psi)/cos(psi); 
    /*vbeam amp contains cos(psi). So replace with sin(psi)
    to get get the gradient */
    
    double delta_psi = delta_v / grad; /*linear extrapolation*/

    
    /* If the cavity is detuned positively, the psi needs to
    be increased to reduce the voltage. Likewise, if the cavity
    is detuned negatively, the psi needs to be decreased to reduce
    the voltage.
    */
        
    int sg = (psi<0) - (psi>0);

    /* This is to avoid setting a value if grad is 0, as then
    delta_psi is inf, which even when multiplied by 0 gives nan
    */
    if (grad!=0.0){
        vgen[2] += sg*delta_psi*phasegain;
    }
}

static void compute_buffer_mean(double *out_array, double *buffer, long windowlength, long buffersize, long numcolumns){

    int c,p,offset;
    offset = buffersize - windowlength;

    for (p=0; p<numcolumns; p++) {
        out_array[p] = 0.0;
    }
    
    for (c=offset; c<buffersize; c++) {
        for (p=0; p<numcolumns; p++) {
            out_array[p] += buffer[2*c+p];
        }
    }
    
    for (p=0; p<numcolumns; p++) {
        out_array[p] /= windowlength ; 
    }
}

int check_buffer_length(double *buffer, long buffersize, long numcolumns){
    int c;
    int bufferlengthnow=0;
    for (c=0; c<numcolumns*buffersize; c++){
        if (buffer[c]!=0.0){
            bufferlengthnow += 1;
        }
    }
    bufferlengthnow /= numcolumns;
    return bufferlengthnow;
}


static void update_vbeam_set(long fbmode, double *vbeam_set,
                             double *vbeamk, double *vbeam_buffer,
                             long buffersize, long windowlength){
    int bufferlengthnow = 0;
    // If FBMode is set to ONETURN, then set the vbeam and move on
    if(fbmode==1){
        vbeam_set[0] = vbeamk[0];
        vbeam_set[1] = vbeamk[1];        
    }
    // If FBMode is set to WINDOW, compute the vbeam_set from the buffer
    
    else if(fbmode==2){
        // Compute the length of the buffer as we will not act until 
        // the buffer is full. (2 arrays of vbeam and psi)
        
        bufferlengthnow = check_buffer_length(vbeam_buffer, buffersize, 2);

        if(bufferlengthnow >= windowlength){
            compute_buffer_mean(vbeam_set, vbeam_buffer, windowlength, buffersize, 2);
        } 
    }
}


