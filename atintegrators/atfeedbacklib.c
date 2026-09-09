#include "atconstants.h"
#include "atelem.c"
#include <math.h>
#include <float.h>
#include <complex.h>
#ifdef MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif



void roll_array(double *arr, int arr_len){


    if(arr_len>1){
        int idx = 0;
        double tmp=0.0;
        double tmp2=0.0;    
        tmp = arr[arr_len-1];
       
        for(idx=0;idx<arr_len;idx++){
            tmp2 = arr[idx];
            arr[idx] = tmp;
            tmp = tmp2;
        }
    }    
}


static void update_vbeam_set(long fbmode, double *vbeam_set,
                             double *vbeamk, double *vbeam_buffer,
                             long buffersize, long windowlength){

    vbeam_set[0] = vbeamk[0];
    vbeam_set[1] = vbeamk[1];        
}


void compute_buffer_mean(double *out_array, double *buffer, long windowlength, long buffersize, long numcolumns){

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



static void compute_set_params(double *vbeam, double *vgen, double phis, double *vgen_set){

    double vbeamr_meas = vbeam[0]*cos(vbeam[1]);
    double vbeami_meas = vbeam[0]*sin(vbeam[1]);
    
    double vgenr_meas = -vgen[0]*sin(vgen[1]);
    double vgeni_meas = vgen[0]*cos(vgen[1]);      
    
    double vcavr_meas = vgenr_meas + vbeamr_meas;
    double vcavi_meas = vgeni_meas + vbeami_meas;   

    double vcav_meas = sqrt(vcavr_meas*vcavr_meas + vcavi_meas*vcavi_meas); 
    double phis_meas = -atan2(vcavr_meas, vcavi_meas);

    double meas_psi = vgen[1] - phis_meas;
    
    if(meas_psi<-TWOPI/2){
        meas_psi += TWOPI;
    }else if(meas_psi > TWOPI/2){
        meas_psi -= TWOPI;
    }
    
    vgen_set[0] = vcav_meas;
    vgen_set[1] = phis_meas;
    vgen_set[2] = meas_psi;

}
static void update_vgen(double *vcav, double *vgen, double *vcav_meas, double voltgain,
                        double phasegain, double *VoltDelay, double *PhaseDelay, int delay){

    double diff_Amp = VoltDelay[delay-1] - vcav[0];
    double diff_Phase = PhaseDelay[delay-1] - vcav[1];
    vgen[0] -= voltgain * diff_Amp;
    vgen[1] -= phasegain * diff_Phase;
    
    roll_array(VoltDelay, delay);
    roll_array(PhaseDelay, delay);
    
    VoltDelay[0] = vcav_meas[0];
    PhaseDelay[0] = vcav_meas[1];    
}

static void compute_tuner(double *vcav_meas, double *vgen_arr,
                          double *TunerParams, double TunerGain, double TunerAveragingPeriod,
                          double TunerOffset){

    if(TunerGain>0){
        TunerParams[0] += 1; // TunerCount        
        TunerParams[1] += (vcav_meas[2] - vgen_arr[2]); //TunerDiff
        
        if(TunerParams[0]==TunerAveragingPeriod){
            TunerParams[1] = (TunerParams[1]/TunerAveragingPeriod) + TunerOffset;
            vgen_arr[2] += TunerGain * TunerParams[1];
            TunerParams[0] = 0.0; //TunerCount
            TunerParams[1] = 0.0; //TunerDiff
        }
    }     
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
