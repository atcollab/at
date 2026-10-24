#include "atconstants.h"
#include "atelem.c"
#include <math.h>
#include <float.h>
#include <complex.h>
#ifdef MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif


void init_IIR(double cutoff, double *IIRcoef, double *IIRout, double T1, int every, double Vc){
    
    if(cutoff==0){
        IIRcoef[0] =  1.0;
    }else{
        double omega = TWOPI * cutoff;
        double T = T1 * every;
        double alpha = cos(omega * T) - 1;
        double tmp = alpha*alpha - 2*alpha;
        if(tmp > 0){
            IIRcoef[0] = alpha + sqrt(tmp);
        }else{
            IIRcoef[0] = T * cutoff * TWOPI;
        }
    }
    IIRout[0] = Vc;
    IIRout[1] = 0.0;
}

void IIR(double complex input, double *IIRcoef, double *IIRout){
    IIRout[0] = (1 - IIRcoef[0]) * IIRout[0] + IIRcoef[0] * creal(input);
    IIRout[1] = (1 - IIRcoef[0]) * IIRout[1] + IIRcoef[0] * cimag(input);
}

double complex Vg2Ig(double complex vgen_phasor, double psi, double RL){
    /*
    Return Ig from Vg (assuming constant Vg).

    Eq.25 of ref [2] assuming the dVg/dt = 0.
    */
    return (vgen_phasor / RL) * (1 - _Complex_I * tan(psi));
}

void init_PIL(double vgen, double thetag, double psi, double RL, double *vbeam,
              double Vc, bool FF, double cutoff, double T1, int every,
              int ring_harmn, int samplenum, int record_size,
              double *ig_phasor_real, double *ig_phasor_imag,
              double *generator_phasor_record_real, double *generator_phasor_record_imag,
              double *vc_previous_real, double *vc_previous_imag,
              double *diff_record_real, double *diff_record_imag,
              double *I_record, double *FFconst, double *IIRcoef, double *IIRout,
              long *samplestart){

    double Vg_real = -vgen*sin(thetag);
    double Vg_imag = vgen*cos(thetag);
    double complex ig = Vg2Ig(Vg_real + _Complex_I * Vg_imag, psi, RL);
    double complex vcav = Vg_real + _Complex_I * Vg_imag + vbeam[0]*cexp(_Complex_I*vbeam[1]);
    int idx;
    for(idx=0;idx<ring_harmn;idx++){
        ig_phasor_real[idx] = creal(ig);
        ig_phasor_imag[idx] = cimag(ig);
        generator_phasor_record_real[idx] = Vg_real;
        generator_phasor_record_imag[idx] = Vg_imag;
    }
    for(idx=0;idx<samplenum;idx++){
        vc_previous_real[idx] = creal(vcav);
        vc_previous_imag[idx] = cimag(vcav);
    }
    for(idx=0;idx<record_size;idx++){
        diff_record_real[idx] = 0.0;
        diff_record_imag[idx] = 0.0;
    }
    FFconst[0] = FF ? creal(ig) : 0.0;
    FFconst[1] = FF ? cimag(ig) : 0.0;
    I_record[0] = 0.0;
    I_record[1] = 0.0;
    samplestart[0] = 0;
    init_IIR(cutoff, IIRcoef, IIRout, T1, every, Vc);
}

void set_cavity_record_phasor_array(double *vbunch, 
                                    double *cavity_phasor_record_real, double *cavity_phasor_record_imag,
                                    double *generator_phasor_record_real, double *generator_phasor_record_imag,
                                    int ring_harmn){
        int idx;
        double complex tmp = 0 + _Complex_I*0;
        
        for(idx=0;idx<ring_harmn;idx++){
            tmp = vbunch[idx]*cexp(_Complex_I*vbunch[idx+ring_harmn]);
            cavity_phasor_record_real[idx] = creal(tmp) + generator_phasor_record_real[idx];
            cavity_phasor_record_imag[idx] = cimag(tmp) + generator_phasor_record_imag[idx];
        }
}

static double complex vc_sample(double *vc_previous_real, double *vc_previous_imag,
                                double *cavity_phasor_record_real, double *cavity_phasor_record_imag,
                                int samplenum, int ibucket){
    if(ibucket < samplenum){
        return vc_previous_real[ibucket] + _Complex_I * vc_previous_imag[ibucket];
    }else{
        return cavity_phasor_record_real[ibucket-samplenum] + _Complex_I * cavity_phasor_record_imag[ibucket-samplenum];
    }
}


void Ig2Vg(double *generator_phasor_record_real, double *generator_phasor_record_imag,
           double *ig_phasor_real, double *ig_phasor_imag,
           double filling_time, double psi,
           double kloss, double T1, int ring_harmn, double *vgen_arr){

    double complex r = cexp(-T1 / filling_time * (1 - _Complex_I * tan(psi)));
    double complex vg = generator_phasor_record_real[ring_harmn-1] + _Complex_I * generator_phasor_record_imag[ring_harmn-1];
    int idx;
    for(idx=0;idx<ring_harmn;idx++){
        vg = r * vg + kloss * T1 * (ig_phasor_real[idx] + _Complex_I * ig_phasor_imag[idx]);
        generator_phasor_record_real[idx] = creal(vg);
        generator_phasor_record_imag[idx] = cimag(vg);
    }

    double mean_vg = 0.0;
    double mean_real = 0.0;
    double mean_thetag=0.0;
    double mean_imag = 0.0;
    for(idx=0;idx<ring_harmn;idx++){
        mean_real += generator_phasor_record_real[idx] / ring_harmn;
        mean_imag += generator_phasor_record_imag[idx] / ring_harmn;
    }

    mean_vg = sqrt(mean_real*mean_real + 
                   mean_imag*mean_imag);
    mean_thetag = -atan2(mean_real, mean_imag);

    vgen_arr[0] = mean_vg;
    vgen_arr[1] = mean_thetag;
}

void update_vc_previous(double *vc_previous_real, double *vc_previous_imag, int samplenum, int ring_harmn, double *cavity_phasor_record_real, double *cavity_phasor_record_imag){
    int idx=0;
    for(idx=0;idx<samplenum;idx++){
        vc_previous_real[idx] = cavity_phasor_record_real[ring_harmn-samplenum+idx];
        vc_previous_imag[idx] = cavity_phasor_record_imag[ring_harmn-samplenum+idx];
    }    
}


void roll_array(double *arr, int arr_len, int shift){
    
    memmove(arr + shift, arr, (arr_len-shift) * sizeof(*arr));
    memset(arr, 0, shift * sizeof(*arr));
}



void compute_set_params(double *vbeam, double *vgen, double phis, double *vcav_meas){

    double vbeamr_meas = vbeam[0]*cos(vbeam[1]);
    double vbeami_meas = vbeam[0]*sin(vbeam[1]);
    
    double vgenr_meas = -vgen[0]*sin(vgen[1]);
    double vgeni_meas = vgen[0]*cos(vgen[1]);
    
    double vcavr_meas = vgenr_meas + vbeamr_meas;
    double vcavi_meas = vgeni_meas + vbeami_meas;   

    double vcav_amp_meas = sqrt(vcavr_meas*vcavr_meas + vcavi_meas*vcavi_meas); 
    double phis_meas = -atan2(vcavr_meas, vcavi_meas);

    double meas_psi = vgen[1] - phis_meas;
    
    // This part is needed to make sure there is not a 2pi
    // phase difference (sometimes seen with harmonic cavity)

    if(meas_psi<-TWOPI/2){
        meas_psi += TWOPI;
    }else if(meas_psi > TWOPI/2){
        meas_psi -= TWOPI;
    }
    
    vcav_meas[0] = vcav_amp_meas;
    vcav_meas[1] = phis_meas;
    vcav_meas[2] = meas_psi;

}
static void update_vgen(double *vcav, double *vgen, double *vcav_meas, double voltgain,
                        double phasegain, double *VoltDelay, double *PhaseDelay, int delay){
    double diff_Amp = VoltDelay[delay-1] - vcav[0];
    double diff_Phase = PhaseDelay[delay-1] - vcav[1];
    vgen[0] -= voltgain * diff_Amp;
    vgen[1] -= phasegain * diff_Phase;
    
    roll_array(VoltDelay, delay, 1);
    roll_array(PhaseDelay, delay, 1);

    
    VoltDelay[0] = vcav_meas[0];
    PhaseDelay[0] = vcav_meas[1];    
}

void compute_tuner(double *vcav_meas, double *vgen_arr,
                          double *TunerParams, double TunerGain, double TunerAveragingPeriod,
                          double TunerOffset){

    TunerParams[0] += 1; // TunerCount        
    TunerParams[1] += (vcav_meas[2] - vgen_arr[2]); //TunerDiff
    
    if(TunerParams[0]==TunerAveragingPeriod){
        TunerParams[1] = (TunerParams[1]/TunerAveragingPeriod) + TunerOffset;
        vgen_arr[2] += TunerGain * TunerParams[1];
        TunerParams[0] = 0.0; //TunerCount
        TunerParams[1] = 0.0; //TunerDiff
    }
}

void update_passive_frequency(double *vbeam, double *vcav, double *vgen,
                                     double *TunerParams, double TunerGain, double TunerAveragingPeriod){
    /* The cavity voltage is
    V(t) = 2*I0*rs*cos(psi)*exp(i(wt+psi))
    We save the amplitude of vbeam, so the exponent goes to 1.
    Therefore vbeam[0] = 2*I0*rs*cos(psi) which is the cavity voltage.
    */

    double vset = vcav[0]; /* desired vbeam */
    double psi = vgen[2]; /*current psi */
    double vpeak = vbeam[0]; /* Peak amplitude of cavity voltage */
    
    /*vbeam amp contains cos(psi). So replace with sin(psi)
    to get get the gradient */
    double grad = vpeak*sin(psi)/cos(psi); 
    double delta_psi = 0.0;

    /* If the cavity is detuned positively, the psi needs to
    be increased to reduce the voltage. Likewise, if the cavity
    is detuned negatively, the psi needs to be decreased to reduce
    the voltage.
    */
        
    int sg = (psi<0) - (psi>0);

    TunerParams[0] += 1; // TunerCount        
    TunerParams[1] += vset - vpeak; //DeltaV
    if(TunerParams[0]==TunerAveragingPeriod){
        TunerParams[1] = (TunerParams[1]/TunerAveragingPeriod);

        delta_psi = TunerParams[1] / grad; /*linear extrapolation*/

        /* This is to avoid setting a value if grad is 0, as then
        delta_psi is inf, which even when multiplied by 0 gives nan
        */
        if (grad!=0.0){
            vgen[2] += sg*delta_psi*TunerGain;
        }
        
        TunerParams[0] = 0.0; //TunerCount
        TunerParams[1] = 0.0; //TunerDiff
    }   
}


void fill_ig_phasor(double *ig_phasor_real, double *ig_phasor_imag,
                    int ring_harmn){
                    
     double fill_value_real = ig_phasor_real[ring_harmn-1];
     double fill_value_imag = ig_phasor_imag[ring_harmn-1];
     int i=0;
     for(i=0;i<ring_harmn;i++){
        ig_phasor_real[i] = fill_value_real;
        ig_phasor_imag[i] = fill_value_imag;
     }                         
                         
}

void track_PIL(double *vc_previous_real, double *vc_previous_imag,
               double *cavity_phasor_record_real, double *cavity_phasor_record_imag,
               double *ig_phasor_real, double *ig_phasor_imag,
               long *samplestart, int samplenum, int record_size,
               double *diff_record_real, double *diff_record_imag,
               double *FFconst, double *gain, double *I_record,
               double Vc, double theta,
               double *generator_phasor_record_real, double *generator_phasor_record_imag,
               double filling_time, double kloss, double T1, int ring_harmn, double *vgen_arr,
               double *IIRout, double *IIRcoef,
               int every,
               double psi, double rshunt,
               int open
               ){
    int index=0;
    int index2=0;

    double diff_real, diff_imag;
    int ibucket;
    int window_start = samplestart[0];
    double complex vc_sum = 0.0;
    double complex mean_vc = 0.0 + _Complex_I * 0.0;
    double complex fb_value, ig;

    for(ibucket=window_start;ibucket<window_start+samplenum;ibucket++){
        vc_sum += vc_sample(vc_previous_real, vc_previous_imag,
                            cavity_phasor_record_real, cavity_phasor_record_imag,
                            samplenum, ibucket);
    }

    fill_ig_phasor(ig_phasor_real, ig_phasor_imag, ring_harmn);
    
    for(index=samplestart[0];index<ring_harmn;index+=every){

        diff_real = diff_record_real[record_size-1];
        diff_imag = diff_record_imag[record_size-1];

        I_record[0] += diff_real * T1;
        I_record[1] += diff_imag * T1;

        fb_value = gain[0] * (diff_real + _Complex_I * diff_imag) 
                 + gain[1] * (I_record[0] + _Complex_I * I_record[1]);
        
        ig = Vg2Ig(fb_value, psi, rshunt);

        for(index2=index;index2<index+every && index2<ring_harmn;index2++){
            ig_phasor_real[index2] = creal(ig) + FFconst[0];
            ig_phasor_imag[index2] = cimag(ig) + FFconst[1];
            }


        roll_array(diff_record_real, record_size, 1);
        roll_array(diff_record_imag, record_size, 1);

        
        for(ibucket=window_start;ibucket<index;ibucket++){
            vc_sum -= vc_sample(vc_previous_real, vc_previous_imag,
                                cavity_phasor_record_real, cavity_phasor_record_imag,
                                samplenum, ibucket);
            vc_sum += vc_sample(vc_previous_real, vc_previous_imag,
                                cavity_phasor_record_real, cavity_phasor_record_imag,
                                samplenum, ibucket+samplenum);
        }
        window_start = index;
        mean_vc = vc_sum / samplenum * cexp(-_Complex_I * (theta+TWOPI/4));
        IIR(mean_vc, IIRcoef, IIRout);

        diff_record_real[0] = Vc - IIRout[0];
        diff_record_imag[0] = -IIRout[1];
    }

    samplestart[0] = index - ring_harmn;

    update_vc_previous(vc_previous_real, vc_previous_imag,
                       samplenum, ring_harmn,
                       cavity_phasor_record_real, cavity_phasor_record_imag);

    if(open==0){
        Ig2Vg(generator_phasor_record_real, generator_phasor_record_imag,
              ig_phasor_real, ig_phasor_imag,
              filling_time, psi,
              kloss, T1, ring_harmn, vgen_arr);
    }
}


