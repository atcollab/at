#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "atimplib.c"
#include "attrackfunc.c"

#ifndef _MSC_VER  
#include "atfeedbacklib.c"
#include <complex.h>
#endif

/*
 * BeamLoadingCavity pass method by Simon White.  
 *
 */

struct elem
{
  int nslice; int nturnsw;
  int cavitymode; int fbmode;
  int buffersize;
  int openloop;
  double normfact;
  double *turnhistory;
  double *z_cuts;
  int delay;  double *VoltDelay; double *PhaseDelay;
  double *gain;
  double TunerOffset; int TunerAveragingPeriod; double TunerGain; double *TunerParams;  
  double Length;
  double Voltage;
  double Energy;
  double Frequency;
  double TimeLag;
  double Qfactor;
  double Rshunt;
  double Beta;
  double phis;
  double ts;
  double *vbunch;
  double *vbeam_phasor; double *vbeam;
  double *vcav; double *vgen;
  double *vgen_buffer; double *vbeam_buffer; double *vbunch_buffer;
  int samplenum; int every; int recordsize;
  int ff; double cutoff;
  double *ig_phasor;
  double *generator_phasor_record; double *cavity_phasor_record;
  double *vc_previous;
  double *diff_record;
  long *samplestart;
  double *I_record; double *FFconst; double *IIRout; double *IIRcoef;
  double *pil_init;
  long dampermode;
  long M0_Delay;
  double *M0_PID_Params;
  double *M0_Gain;
  double *M0_Buffer;
  double M0_VRatio;
  }; 


void write_buffer(double *data, double *buffer, int datasize, int buffersize){
    if(buffersize>1){
        memmove(buffer, buffer + datasize, datasize*(buffersize-1)*sizeof(double));
    }
    memcpy(buffer + datasize*(buffersize-1), data, datasize*sizeof(double));
}
   

void compute_mean_dp(double *r_in, int num_particles, double *mdp, double *npart){
    int c = 0;
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        double *r6 = r_in+c*6;

        if (!atIsNaN(r6[0])) {
            npart[0] += 1; 
            mdp[0] += r6[4]; //dp
        }
    }

    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE, &npart, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &mdp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
    
    if (npart[0]>0.0){
        mdp[0] /= npart[0];
    }

}


void BeamLoadingCavityPass(double *r_in, int num_particles, int nbunch,
                           double *bunch_spos, double *bunch_currents, 
                           double *fillpattern,
                           double circumference,
                           int nturn, double energy, int harmonic_number,
                           int iturn,
                           struct elem *Elem) {
  
    // some wake definitions
    long nslice = Elem->nslice;
    long nturnsw = Elem->nturnsw; 
    long buffersize = Elem->buffersize;
    double normfact = Elem->normfact;  
    
    // load buffer pointers
    double *turnhistory = Elem->turnhistory;
    double *vgen_buffer = Elem->vgen_buffer;
    double *vbeam_buffer = Elem->vbeam_buffer;
    double *vbunch_buffer = Elem->vbunch_buffer;

    double *z_cuts = Elem->z_cuts;
    double *vbunch = Elem->vbunch;
    double *vbeam = Elem->vbeam;
    double *vcav_set = Elem->vcav; /* Vcav set points amplitude, phase */
    double *vbeam_phasor = Elem->vbeam_phasor;
    double *vgen_arr = Elem->vgen; // [vgen, thetag, psi]
        
    // ring parameters
    int ring_harmn = harmonic_number; // Ring harmonic number (number of buckets)
    
    // cavity related parameters
    double le = Elem->Length;
    double rffreq = Elem->Frequency;
    int harmn = rffreq * circumference / C0 ;    // cavity harmonic number

    double tlag = Elem->TimeLag;
    double qfactor = Elem->Qfactor; //loaded
    double rshunt = Elem->Rshunt; // loaded 
    double beta = Elem->Beta; //relativist beta not cavity beta
    double ts = Elem->ts; // expected beam final position
    double f0 = beta *C0 / circumference; // revolution frequency
    
    long cavitymode = Elem->cavitymode;
    long fbmode = Elem->fbmode;


    //Tuner Variables
    double TunerGain = Elem->TunerGain;
    double TunerOffset = Elem->TunerOffset;
    int TunerAveragingPeriod = Elem->TunerAveragingPeriod;
    double *TunerParams = Elem->TunerParams; //TunerParams[0] is TunerCount, TunerParams[1] is TunerDiff

    //if fb mode is PROP then gain[0] is Voltgain and gain[1] is PhaseGain
    //if fb mode is PROP_INTEGRAL then gain[0] is Prop gain and gain[1] is integral gain
    double *gain = Elem->gain;

    //if fb mode is PROP then delay is in units of turns
    //if fb mode is PROP_INTEGRAL then delay is in units of buckets    
    int delay = Elem->delay; 
    
    // Only used for fb mode PROP
    double *VoltDelay = Elem->VoltDelay;
    double *PhaseDelay = Elem->PhaseDelay;


            
    double cutoff = Elem->cutoff;  // cutoff frequency 

    int samplenum = Elem->samplenum; // from 0  to samplenum buckets...
    int every = Elem->every;  // ...in steps of every 
    int FF = Elem->ff; //Use the feedforward constant? 
    int record_size = Elem->recordsize; // overlap coming from delay and every sampling
    int open = Elem->openloop; // do you want to apply the correction?
    
    double *I_record = Elem->I_record; // real and imaginary of the integral part of the loop            
    double *FFconst = Elem->FFconst; // the feedfoward constant to use
    double *IIRout = Elem->IIRout;
    double *IIRcoef = Elem->IIRcoef;
    double *pil_init = Elem->pil_init;
    
    /* Here we have to declare pointers for the PI Loop
    They have to be defined outside of an if statement*/
    
    double *ig_phasor_real = Elem->ig_phasor;
    double *ig_phasor_imag = Elem->ig_phasor + ring_harmn;
                      
    double *generator_phasor_record_real = Elem->generator_phasor_record;
    double *generator_phasor_record_imag = Elem->generator_phasor_record + ring_harmn;
    double *cavity_phasor_record_real = Elem->cavity_phasor_record;
    double *cavity_phasor_record_imag = Elem->cavity_phasor_record + ring_harmn;
    double *vc_previous_real = Elem->vc_previous; 
    double *vc_previous_imag = Elem->vc_previous + samplenum;
    double *diff_record_real = Elem->diff_record; 
    double *diff_record_imag = Elem->diff_record + record_size;
    long *samplestart = Elem->samplestart;

    // End of loading element attributes.

    // Load M0 Damper Attributes
    long dampermode = Elem->dampermode;
    double *M0_Gain = Elem->M0_Gain;
    double *M0_PID_Params = Elem->M0_PID_Params;
    long M0_Delay = Elem->M0_Delay;
    double *M0_Buffer = Elem->M0_Buffer;
    double M0_VRatio = Elem->M0_VRatio;
    double mdp[] = {0.0};
    double npart[] = {0.0};
    //
    
    double vbeam_set[] = {vbeam[0], vbeam[1]}; 
    double vcav_meas[] = {0.0, 0.0, 0.0};

    double tot_current = 0.0;

    
    int i, c;
    size_t sz = nslice*nbunch*sizeof(double) + num_particles*sizeof(int);
    
    // Two empty pointers that will be pointed to later.
    int *pslice;
    double *vbeam_kicks; /* This used to be kz, it is the kick that is applied */

    // For now it is easier to use the variable, can be removed later
    double vgen = vgen_arr[0];
    double gen_phase = vgen_arr[1];
    double psi = vgen_arr[2];
    double delta = pow(rffreq * tan(psi) / qfactor, 2) + 4 * pow(rffreq,2);
    double freqres = (rffreq * tan(psi) / qfactor + sqrt(delta)) / 2;

    double tot_lag_phase = (tlag+ts)*rffreq*TWOPI/C0;

    // parameters needed to computed cavity response
    double filling_time = 2*qfactor / (TWOPI * freqres);
    double T1 = circumference / (ring_harmn * beta * C0);
    double kloss = rshunt * TWOPI * freqres / (2 * qfactor);


    
    for(i=0;i<nbunch;i++){
        tot_current += bunch_currents[i];
    }
    
    double m0_lagshift = 0.0;    
    #ifndef _MSC_VER
    // Mode 0 computation
    if(dampermode==2){
        double error_before = M0_PID_Params[1];
        double total_voltage = vcav_set[0] / M0_VRatio;
        double M0_U0 = total_voltage * sin(vcav_set[1]);
        compute_mean_dp(r_in, num_particles, mdp, npart);
        double dp_eV = mdp[0] * energy;
        double new_phis = -asin((M0_U0 - dp_eV)/total_voltage);
        double new_timelag = -new_phis*C0/(TWOPI*rffreq);        
        double delta = new_timelag - tlag;
        M0_PID_Params[0] = delta;
        M0_PID_Params[1] += delta/f0;
    
        double derivative = (M0_PID_Params[1] - error_before)*f0;
    
        double fb_value = M0_PID_Params[0]*M0_Gain[0] + M0_PID_Params[1]*M0_Gain[1] + M0_PID_Params[2]*M0_Gain[2];
        roll_array(M0_Buffer, M0_Delay, 1);
        M0_Buffer[M0_Delay-1] = fb_value;
        m0_lagshift = M0_Buffer[0];
    }
    #endif
    //
    
    
    /*Track RF cavity is always done. */
    trackRFCavity(r_in, le, vgen/energy, rffreq, harmn, tlag + m0_lagshift, -gen_phase - tot_lag_phase, nturn, circumference/C0, num_particles);

    #ifndef _MSC_VER
    /*Only allocate memory if current is > 0*/
    if(tot_current>0 && rshunt > 0){
        // allocate memory and set pointers
        void *buffer = atMalloc(sz);        
        double *dptr = (double *) buffer;
        int *iptr;
        vbeam_kicks = dptr;
        dptr += nslice*nbunch;
        iptr = (int *) dptr;
        pslice = iptr; 
        iptr += num_particles;


        rotate_table_history(nturnsw, nslice*nbunch, turnhistory, circumference);
        slice_bunch(r_in, num_particles, nslice, nturnsw, nbunch, bunch_spos,
                    bunch_currents, turnhistory, pslice, z_cuts);
        compute_kicks_phasor(nslice, nbunch, nturnsw, turnhistory, normfact, vbeam_kicks,
                             freqres, qfactor, rshunt, vbeam_phasor, circumference, energy,
                             beta, vbeam_set, vbunch, bunch_spos, ring_harmn, fillpattern, ts);             

        /*apply kicks*/
        for (c=0; c<num_particles; c++) {
            double *r6 = r_in+c*6;
            int islice=pslice[c];
            if (!atIsNaN(r6[0])) {         
                r6[4] += vbeam_kicks[islice];                 
            }
        }
        
        // the beam tracking is now done, from now on we compute the feedbacks
               
        // First write the values to the buffer
        if(buffersize>0){
            write_buffer(vbeam_set, vbeam_buffer, 2, buffersize);
            write_buffer(vgen_arr, vgen_buffer, 3, buffersize);
            write_buffer(vbunch, vbunch_buffer, 2*ring_harmn, buffersize);
        }   


        compute_set_params(vbeam_set, vgen_arr, vcav_set[1], vcav_meas);
        
        if(cavitymode==1){
            // If CavityMode=ACTIVE
            if(fbmode==1){
                // If FBMode=PROP
                update_vgen(vcav_set, vgen_arr, vcav_meas, gain[0], gain[1], VoltDelay, PhaseDelay, delay, open);
            }
            if(fbmode==2){
                // If FBMode=PROP_INTEGRAL
                if(pil_init[0]==0){
                    init_PIL(vgen, gen_phase, psi, rshunt, vbeam_set,
                             vcav_set[0], FF, cutoff, T1, every,
                             ring_harmn, samplenum, record_size,
                             ig_phasor_real, ig_phasor_imag,
                             generator_phasor_record_real, generator_phasor_record_imag,
                             vc_previous_real, vc_previous_imag,
                             diff_record_real, diff_record_imag,
                             I_record, FFconst, IIRcoef, IIRout, samplestart);
                    pil_init[0] = 1.0;
                }

                set_cavity_record_phasor_array(vbunch,
                                               cavity_phasor_record_real, cavity_phasor_record_imag,
                                               generator_phasor_record_real, generator_phasor_record_imag,
                                               ring_harmn); 

                track_PIL(vc_previous_real, vc_previous_imag,
                          cavity_phasor_record_real, cavity_phasor_record_imag,
                          ig_phasor_real, ig_phasor_imag,
                          samplestart, samplenum, record_size,
                          diff_record_real, diff_record_imag,
                          FFconst, gain, I_record,
                          vcav_set[0], vcav_set[1],
                          generator_phasor_record_real, generator_phasor_record_imag,
                          filling_time, kloss, T1, ring_harmn, vgen_arr,
                          IIRout, IIRcoef,
                          every,
                          psi, rshunt,
                          open
                          );    

            }
            
        }else if(cavitymode==3){     
            /// If CavityMode=PASSIVE_VOLTAGE
            update_passive_frequency(vbeam_set, vcav_set, vgen_arr,
                                     TunerParams, TunerGain, TunerAveragingPeriod);
        }


        // Here is where the tuner is calculated and applied for PROP and PROP_INTEGRAL 
        // If TunerGain is zero, it is skipped 
        // If you are in PASSIVE_VOLTAGE, then the psi 
        //  has already been updated so don't need to come here again. 
        if(TunerGain>0 && cavitymode==1){
            compute_tuner(vcav_meas, vgen_arr,
                          TunerParams, TunerGain, TunerAveragingPeriod, TunerOffset);
        }            

        vbeam[0] = vbeam_set[0];
        vbeam[1] = vbeam_set[1];
        atFree(buffer);
    }
    #endif
}


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double rl = Param->RingLength;
    double energy;
    int nturn = Param->nturn;
    if (!Elem) {
        long nslice, nturns, cavitymode, fbmode, dampermode, buffersize;
        long delay, every, samplenum, ff, recordsize, openloop;
        double TunerGain, TunerOffset, TunerAveragingPeriod, *TunerParams;
        double *VoltDelay, *PhaseDelay, *pil_init;
        
        double wakefact, Energy, Frequency, TimeLag, Length;
        double normfact, qfactor, rshunt, beta, phis, ts, cutoff;
        double *gain;
        double *turnhistory;
        double *vgen_buffer;
        double *vbeam_buffer;
        double *vbunch_buffer;
        double *z_cuts;

        double *vbunch;
        double *vbeam_phasor;
        double *vbeam;
        double *vgen;
        double *vcav;

        double *ig_phasor;
        double *generator_phasor_record;
        double *cavity_phasor_record;
        double *vc_previous;
        double *diff_record;
        long *samplestart;
        double *I_record;
        double *FFconst;
        double *IIRcoef;
        double *IIRout;
        
        double *M0_Gain, *M0_Buffer, *M0_PID_Params;
        long M0_Delay;
        double M0_VRatio;
        
        /*attributes for RF cavity*/
        Length=atGetDouble(ElemData,"Length"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        
        /*attributes for resonator*/
        nslice=atGetLong(ElemData,"_nslice"); check_error();
        nturns=atGetLong(ElemData,"_nturns"); check_error();
        buffersize=atGetLong(ElemData,"_buffersize"); check_error();
        cavitymode=atGetLong(ElemData,"_cavitymode"); check_error();
        fbmode=atGetLong(ElemData,"_fbmode"); check_error();
        wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
        qfactor=atGetDouble(ElemData,"Qfactor"); check_error();
        rshunt=atGetDouble(ElemData,"Rshunt"); check_error();
        beta=atGetDouble(ElemData,"_beta"); check_error();
        normfact=atGetDouble(ElemData,"NormFact"); check_error();
        gain=atGetDoubleArray(ElemData,"Gain"); check_error();        
        TunerGain=atGetDouble(ElemData,"TunerGain"); check_error();
        turnhistory=atGetDoubleArray(ElemData,"_turnhistory"); check_error();
        vbunch=atGetDoubleArray(ElemData,"_vbunch"); check_error();
        vbeam=atGetDoubleArray(ElemData,"_vbeam"); check_error();
        vcav=atGetDoubleArray(ElemData,"_vcav"); check_error();
        vgen=atGetDoubleArray(ElemData,"_vgen"); check_error();
        vbeam_phasor=atGetDoubleArray(ElemData,"_vbeam_phasor"); check_error(); 
        vgen_buffer=atGetDoubleArray(ElemData,"_vgen_buffer"); check_error();
        vbeam_buffer=atGetDoubleArray(ElemData,"_vbeam_buffer"); check_error();
        vbunch_buffer=atGetDoubleArray(ElemData,"_vbunch_buffer"); check_error();
        
        phis=atGetDouble(ElemData,"_phis"); check_error();
        ts=atGetDouble(ElemData,"_ts"); check_error();
        
        openloop=atGetLong(ElemData,"OpenLoop"); check_error();
        dampermode=atGetLong(ElemData,"_dampermode"); check_error();
        /*optional attributes*/        
        delay=atGetOptionalLong(ElemData,"_delay",1); check_error();
        
        VoltDelay=atGetOptionalDoubleArray(ElemData,"_voltdelay"); check_error();
        PhaseDelay=atGetOptionalDoubleArray(ElemData,"_phasedelay"); check_error();

        every=atGetOptionalLong(ElemData,"_every",1); check_error();
        samplenum=atGetOptionalLong(ElemData,"_samplenum",1); check_error();        
        cutoff=atGetOptionalDouble(ElemData,"cutoff",0); check_error();        
        ff=atGetOptionalLong(ElemData,"FF",1); check_error();
        recordsize=atGetOptionalLong(ElemData,"_recordsize",1); check_error();

        ig_phasor=atGetOptionalDoubleArray(ElemData,"_ig_phasor"); check_error();
        generator_phasor_record=atGetOptionalDoubleArray(ElemData,"_generator_phasor_record"); check_error();
        cavity_phasor_record=atGetOptionalDoubleArray(ElemData,"_cavity_phasor_record"); check_error();

        vc_previous=atGetOptionalDoubleArray(ElemData,"_vc_previous"); check_error();
        diff_record=atGetOptionalDoubleArray(ElemData,"_diff_record"); check_error();        
        samplestart=atGetOptionalLongArray(ElemData,"_samplestart"); check_error();        
        I_record=atGetOptionalDoubleArray(ElemData,"_I_record"); check_error();
        FFconst=atGetOptionalDoubleArray(ElemData,"_FFconst"); check_error();
        IIRcoef=atGetOptionalDoubleArray(ElemData,"_IIRcoef"); check_error();
        IIRout=atGetOptionalDoubleArray(ElemData,"_IIRout"); check_error();
        pil_init=atGetOptionalDoubleArray(ElemData,"_pil_init"); check_error();
        
                
        Energy=atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
        
        TunerOffset=atGetOptionalDouble(ElemData,"TunerOffset", 0.0); check_error();
        TunerAveragingPeriod=atGetOptionalDouble(ElemData,"TunerAveragingPeriod",1); check_error();
        TunerParams=atGetOptionalDoubleArray(ElemData,"_TunerParams"); check_error();
        

        M0_Gain=atGetOptionalDoubleArray(ElemData,"M0_Gain"); check_error();
        M0_Delay=atGetOptionalLong(ElemData,"_M0_Delay",1); check_error();
        M0_Buffer=atGetOptionalDoubleArray(ElemData,"_M0_Buffer"); check_error();
        M0_PID_Params=atGetOptionalDoubleArray(ElemData,"_M0_PID_Params"); check_error();
        M0_VRatio=atGetOptionalDouble(ElemData,"M0_VRatio", 1.0); check_error();
        
        /* Check energy */
        Energy = atEnergy(Param->energy, Energy); check_error();

        int dimsth[] = {Param->nbunch*nslice*nturns, 4};
        atCheckArrayDims(ElemData,"_turnhistory", 2, dimsth); check_error();
        int dimsvb[] = {Param->harmonic_number, 2};
        atCheckArrayDims(ElemData,"_vbunch", 2, dimsvb); check_error();
       
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        
        Elem->Length=Length;
        Elem->Frequency=Frequency;
        Elem->Energy = Energy;
        Elem->TimeLag=TimeLag;   
        Elem->nslice=nslice;
        Elem->nturnsw=nturns;
        Elem->normfact=normfact*wakefact;
        Elem->turnhistory=turnhistory;
        Elem->Qfactor = qfactor;
        Elem->Rshunt = rshunt;
        Elem->Beta = beta;
        Elem->z_cuts = z_cuts;
        Elem->vbunch = vbunch;
        Elem->vbeam = vbeam;
        Elem->vgen = vgen;
        Elem->vcav = vcav;

        Elem->TunerGain = TunerGain;
        Elem->TunerOffset = TunerOffset;
        Elem->TunerAveragingPeriod = TunerAveragingPeriod;
        Elem->TunerParams = TunerParams;
        Elem->openloop = openloop;
        Elem->gain = gain;
        Elem->delay=delay;
        Elem->VoltDelay=VoltDelay;
        Elem->PhaseDelay=PhaseDelay;  

        Elem->vbeam_phasor = vbeam_phasor;
        Elem->cavitymode = cavitymode;
        Elem->buffersize = buffersize;
        Elem->vgen_buffer = vgen_buffer;
        Elem->vbeam_buffer = vbeam_buffer;
        Elem->vbunch_buffer = vbunch_buffer;

        Elem->fbmode = fbmode;
        Elem->phis = phis;
        Elem->ts = ts;

        Elem->every=every;
        Elem->samplenum=samplenum;
        Elem->cutoff=cutoff;
        Elem->ff=ff;
        Elem->recordsize=recordsize;
        Elem->I_record=I_record;
        Elem->ig_phasor=ig_phasor;
        Elem->generator_phasor_record=generator_phasor_record;
        Elem->cavity_phasor_record=cavity_phasor_record;
        Elem->vc_previous=vc_previous;
        Elem->diff_record=diff_record;
        Elem->samplestart=samplestart;
        Elem->FFconst=FFconst;
        Elem->IIRcoef=IIRcoef;
        Elem->IIRout=IIRout;
        Elem->pil_init=pil_init;

        Elem->dampermode=dampermode;
        Elem->M0_Gain = M0_Gain;
        Elem->M0_Delay = M0_Delay;
        Elem->M0_Buffer = M0_Buffer;
        Elem->M0_PID_Params = M0_PID_Params;
        Elem->M0_VRatio = M0_VRatio;
    }
    energy = atEnergy(Param->energy, Elem->Energy); check_error();

    if(num_particles<Param->nbunch){
        atError("Number of particles has to be greater or equal to the number of bunches."); check_error();
    }else if (num_particles%Param->nbunch!=0){
        atWarning("Number of particles not a multiple of the number of bunches: uneven bunch load."); check_error();
    }
    if(Elem->cavitymode==0 || Elem->cavitymode>=4){
        atError("Unknown cavitymode provided."); check_error();
    } 
    if(Elem->fbmode>=3){
        atError("Unknown fbmode provided."); check_error();
    } 
    if(Elem->fbmode>=3){
        atError("Unknown dampermode provided."); check_error();
    } 
    

    #ifdef _MSC_VER
    atError("Beam loading module not implemented in Windows."); check_error();
    #endif
    
    BeamLoadingCavityPass(r_in,num_particles,Param->nbunch,Param->bunch_spos,
                          Param->bunch_currents, Param->fillpattern, rl, 
                          nturn, energy, Param->harmonic_number, Param->nturn, Elem);
    return Elem;
}

MODULE_DEF(BeamLoadingCavityPass)       /* Dummy module initialisation */
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if(nrhs >= 2) {
      double rest_energy = 0.0;
      double charge = -1.0;
      double *r_in;
      const mxArray *ElemData = prhs[0];
      int num_particles = mxGetN(prhs[1]);
      struct elem El, *Elem=&El;
      long nslice, nturns, cavitymode, fbmode, dampermode, buffersize;
      long delay, every, samplenum, ff, recordsize, openloop;

      double TunerGain, TunerOffset, TunerAveragingPeriod, *TunerParams;
      double *VoltDelay, *PhaseDelay, *pil_init;
      double wakefact, Energy, Frequency, TimeLag, Length;
      double normfact, qfactor, rshunt, beta, phis, ts, cutoff;
      double *gain;
      double *turnhistory;
      double *z_cuts;
      double *vbunch;
      double *vbeam_phasor;
      double *vbeam;
      double *vgen;
      double *vcav;
      double *vgen_buffer;
      double *vbeam_buffer;
      double *vbunch_buffer;
      double *I_record;
      double *FFconst;
      double *IIRcoef;
      double *IIRout;
      
      
      double *ig_phasor;
      double *generator_phasor_record;
      double *cavity_phasor_record;
      double *vc_previous;
      double *diff_record;
      long *samplestart;

      double *M0_Gain, *M0_Buffer, *M0_PID_Params;
      long M0_Delay;
      double M0_VRatio;
      
      /*attributes for RF cavity*/
      Length=atGetDouble(ElemData,"Length"); check_error();
      Frequency=atGetDouble(ElemData,"Frequency"); check_error();
      TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
      /*attributes for resonator*/
      nslice=atGetLong(ElemData,"_nslice"); check_error();
      nturns=atGetLong(ElemData,"_nturns"); check_error();
      buffersize=atGetLong(ElemData,"_buffersize"); check_error();
      cavitymode=atGetLong(ElemData,"_cavitymode"); check_error();
      fbmode=atGetLong(ElemData,"_fbmode"); check_error();
      wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
      qfactor=atGetDouble(ElemData,"Qfactor"); check_error();
      rshunt=atGetDouble(ElemData,"Rshunt"); check_error();
      beta=atGetDouble(ElemData,"_beta"); check_error();
      normfact=atGetDouble(ElemData,"NormFact"); check_error();

      gain=atGetDoubleArray(ElemData,"Gain"); check_error();
      turnhistory=atGetDoubleArray(ElemData,"_turnhistory"); check_error();
      vbunch=atGetDoubleArray(ElemData,"_vbunch"); check_error();
      vbeam=atGetDoubleArray(ElemData,"_vbeam"); check_error();
      vcav=atGetDoubleArray(ElemData,"_vcav"); check_error();
      vgen=atGetDoubleArray(ElemData,"_vgen"); check_error();
      vbeam_phasor=atGetDoubleArray(ElemData,"_vbeam_phasor"); check_error(); 
      vgen_buffer=atGetDoubleArray(ElemData,"_vgen_buffer"); check_error();
      vbeam_buffer=atGetDoubleArray(ElemData,"_vbeam_buffer"); check_error();
      vbunch_buffer=atGetDoubleArray(ElemData,"_vbunch_buffer"); check_error();
      phis=atGetDouble(ElemData,"_phis"); check_error();
      ts=atGetDouble(ElemData,"_ts"); check_error();
      openloop=atGetLong(ElemData,"OpenLoop"); check_error();
      dampermode=atGetLong(ElemData,"_dampermode"); check_error();
      /*optional attributes*/
      Energy=atGetOptionalDouble(ElemData,"Energy",0.0); check_error();
      z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();

      TunerGain=atGetDouble(ElemData,"TunerGain"); check_error();
      TunerAveragingPeriod=atGetOptionalLong(ElemData,"TunerAveragingPeriod",1); check_error();
      TunerOffset=atGetOptionalDouble(ElemData,"TunerOffset",0.0); check_error();
      TunerParams=atGetOptionalDoubleArray(ElemData,"_TunerParams"); check_error();
      
      ig_phasor=atGetOptionalDoubleArray(ElemData,"_ig_phasor"); check_error();
      generator_phasor_record=atGetOptionalDoubleArray(ElemData,"_generator_phasor_record"); check_error();
      cavity_phasor_record=atGetOptionalDoubleArray(ElemData,"_cavity_phasor_record"); check_error();


      delay=atGetOptionalLong(ElemData,"_delay",1); check_error();
      VoltDelay=atGetOptionalDoubleArray(ElemData,"_voltdelay"); check_error();
      PhaseDelay=atGetOptionalDoubleArray(ElemData,"_phasedelay"); check_error();  
          
      every=atGetOptionalLong(ElemData,"_every",1); check_error();
      samplenum=atGetOptionalLong(ElemData,"_samplenum",1); check_error();        
      cutoff=atGetOptionalDouble(ElemData,"cutoff",0); check_error();        
      ff=atGetOptionalLong(ElemData,"FF",1); check_error();
      recordsize=atGetOptionalLong(ElemData,"_recordsize",1); check_error();  
        
      vc_previous=atGetOptionalDoubleArray(ElemData,"_vc_previous"); check_error();
      diff_record=atGetOptionalDoubleArray(ElemData,"_diff_record"); check_error();        
      samplestart=atGetOptionalLongArray(ElemData,"_samplestart"); check_error();        
      I_record=atGetOptionalDoubleArray(ElemData,"_I_record"); check_error();
      FFconst=atGetOptionalDoubleArray(ElemData,"_FFconst"); check_error();
      IIRcoef=atGetOptionalDoubleArray(ElemData,"_IIRcoef"); check_error();
      IIRout=atGetOptionalDoubleArray(ElemData,"_IIRout"); check_error();
      pil_init=atGetOptionalDoubleArray(ElemData,"_pil_init"); check_error();


      M0_Gain=atGetOptionalDoubleArray(ElemData,"M0_Gain"); check_error();
      M0_Delay=atGetOptionalLong(ElemData,"_M0_Delay", 1); check_error();
      M0_Buffer=atGetOptionalDoubleArray(ElemData,"_M0_Buffer"); check_error();
      M0_PID_Params=atGetOptionalDoubleArray(ElemData,"_M0_PID_Params"); check_error();
      M0_VRatio=atGetOptionalDouble(ElemData,"M0_VRatio", 1.0); check_error();
      
      Elem = (struct elem*)atMalloc(sizeof(struct elem));
      Elem->Length=Length;
      Elem->cavitymode=cavitymode;
      Elem->Frequency=Frequency;
      Elem->Energy = Energy;
      Elem->TimeLag=TimeLag;   
      Elem->nslice=nslice;
      Elem->nturnsw=nturns;
      Elem->normfact=normfact*wakefact;
      Elem->turnhistory=turnhistory;
      Elem->Qfactor = qfactor;
      Elem->Rshunt = rshunt;
      Elem->Beta = beta;
      Elem->z_cuts=z_cuts;
      Elem->vbunch = vbunch;
      Elem->vbeam = vbeam;
      Elem->vgen = vgen;
      Elem->vcav = vcav;

      Elem->gain = gain;
      Elem->openloop = openloop;
      
      Elem->TunerGain = TunerGain;
      Elem->TunerOffset = TunerOffset;
      Elem->TunerAveragingPeriod = TunerAveragingPeriod;
      Elem->TunerParams = TunerParams;
      
      Elem->VoltDelay=VoltDelay;
      Elem->PhaseDelay=PhaseDelay;
      Elem->delay=delay;
      
      Elem->vbeam_phasor = vbeam_phasor;
      Elem->buffersize = buffersize;
      Elem->vgen_buffer = vgen_buffer;
      Elem->vbeam_buffer = vbeam_buffer;
      Elem->vbunch_buffer = vbunch_buffer;

      Elem->phis = phis;
      Elem->ts = ts;

      Elem->every=every;
      Elem->samplenum=samplenum;
      Elem->cutoff=cutoff;
      Elem->ff=ff;
      Elem->recordsize=recordsize;
      
      Elem->ig_phasor=ig_phasor;
      Elem->generator_phasor_record=generator_phasor_record;
      Elem->cavity_phasor_record=cavity_phasor_record;
      Elem->vc_previous=vc_previous;
      Elem->diff_record=diff_record;
      Elem->samplestart=samplestart;
      Elem->I_record=I_record;
      Elem->FFconst=FFconst;
      Elem->pil_init=pil_init;
      Elem->IIRcoef=IIRcoef;
      Elem->IIRout=IIRout;
      
      Elem->fbmode = fbmode;

      Elem->dampermode = dampermode;
      Elem->M0_Gain = M0_Gain;
      Elem->M0_Delay = M0_Delay; 
      Elem->M0_Buffer = M0_Buffer;
      Elem->M0_PID_Params = M0_PID_Params;
      Elem->M0_VRatio = M0_VRatio;

      if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);

      if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
      /* ALLOCATE memory for the output array of the same size as the input  */
      plhs[0] = mxDuplicateArray(prhs[1]);
      r_in = mxGetDoubles(plhs[0]);

      double bspos = 0.0;
      double bcurr = 0.0;
      double fillp = 0.0;
      BeamLoadingCavityPass(r_in, num_particles, 1, &bspos, &bcurr, &fillp, 1, 0, Energy, 1, 1, Elem);
  }
  else if (nrhs == 0)
  {   /* return list of required fields */

      plhs[0] = mxCreateCellMatrix(28,1);
      mxSetCell(plhs[0],0,mxCreateString("Length"));
      mxSetCell(plhs[0],1,mxCreateString("Energy"));
      mxSetCell(plhs[0],2,mxCreateString("Frequency"));
      mxSetCell(plhs[0],3,mxCreateString("_nslice"));
      mxSetCell(plhs[0],4,mxCreateString("_nturns"));
      mxSetCell(plhs[0],5,mxCreateString("_cavitymode"));
      mxSetCell(plhs[0],6,mxCreateString("_fbmode"));      
      mxSetCell(plhs[0],7,mxCreateString("_wakefact"));
      mxSetCell(plhs[0],8,mxCreateString("Qfactor"));
      mxSetCell(plhs[0],9,mxCreateString("Rshunt"));
      mxSetCell(plhs[0],10,mxCreateString("_beta"));
      mxSetCell(plhs[0],11,mxCreateString("NormFact"));
      mxSetCell(plhs[0],12,mxCreateString("Gain"));
      mxSetCell(plhs[0],13,mxCreateString("_turnhistory"));
      mxSetCell(plhs[0],14,mxCreateString("_vbunch"));
      mxSetCell(plhs[0],15,mxCreateString("_vbeam"));
      mxSetCell(plhs[0],16,mxCreateString("_vcav"));
      mxSetCell(plhs[0],17,mxCreateString("_vgen"));
      mxSetCell(plhs[0],18,mxCreateString("_vbeam_phasor"));
      mxSetCell(plhs[0],19,mxCreateString("_vgen_buffer"));
      mxSetCell(plhs[0],20,mxCreateString("_vbeam_buffer"));
      mxSetCell(plhs[0],21,mxCreateString("_vbunch_buffer"));
      mxSetCell(plhs[0],22,mxCreateString("_buffersize"));
      mxSetCell(plhs[0],23,mxCreateString("_phis"));
      mxSetCell(plhs[0],24,mxCreateString("_ts"));
      mxSetCell(plhs[0],25,mxCreateString("OpenLoop"));
      mxSetCell(plhs[0],26,mxCreateString("TunerGain"));
      mxSetCell(plhs[0],27,mxCreateString("_dampermode"));
                                      
      if(nlhs>1) /* optional fields */
      {
          plhs[1] = mxCreateCellMatrix(29,1);
          mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
          mxSetCell(plhs[1],1,mxCreateString("ZCuts"));
          mxSetCell(plhs[1],2,mxCreateString("TunerOffset"));
          mxSetCell(plhs[1],3,mxCreateString("TunerAveragingPeriod"));
          mxSetCell(plhs[1],4,mxCreateString("_delay"));
          mxSetCell(plhs[1],5,mxCreateString("_every"));
          mxSetCell(plhs[1],6,mxCreateString("_samplenum"));
          mxSetCell(plhs[1],7,mxCreateString("cutoff"));
          mxSetCell(plhs[1],8,mxCreateString("FF"));
          mxSetCell(plhs[1],9,mxCreateString("_ig_phasor"));
          mxSetCell(plhs[1],10,mxCreateString("_generator_phasor_record"));
          mxSetCell(plhs[1],11,mxCreateString("_cavity_phasor_record"));
          mxSetCell(plhs[1],12,mxCreateString("_vc_previous"));
          mxSetCell(plhs[1],13,mxCreateString("_diff_record"));
          mxSetCell(plhs[1],14,mxCreateString("_samplestart"));
          mxSetCell(plhs[1],15,mxCreateString("_I_record"));
          mxSetCell(plhs[1],16,mxCreateString("_FFconst"));
          mxSetCell(plhs[1],17,mxCreateString("_IIRcoef"));
          mxSetCell(plhs[1],18,mxCreateString("_IIRout"));
          mxSetCell(plhs[1],19,mxCreateString("_recordsize"));
          mxSetCell(plhs[1],20,mxCreateString("_TunerParams"));
          mxSetCell(plhs[1],21,mxCreateString("_voltdelay"));
          mxSetCell(plhs[1],22,mxCreateString("_phasedelay"));
          mxSetCell(plhs[1],23,mxCreateString("_pil_init"));
          mxSetCell(plhs[1],24,mxCreateString("M0_Delay"));
          mxSetCell(plhs[1],25,mxCreateString("_M0_Buffer"));
          mxSetCell(plhs[1],26,mxCreateString("_M0_PID_Params"));
          mxSetCell(plhs[1],27,mxCreateString("M0_Gain"));
          mxSetCell(plhs[1],28,mxCreateString("M0_VRatio"));

      }
  }
  else
  {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
  }
  
}
#endif /* MATLAB_MEX_FILE */

