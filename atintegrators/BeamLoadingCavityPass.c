#include "atconstants.h"
#include "atelem.c"
#include "atimplib.c"
#include "atfeedbacklib.c"
#include "attrackfunc.c"
#include <complex.h>

/*
 * BeamLoadingCavity pass method by Simon White.  
 *
 */

struct elem
{
  int nslice;
  int nturnsw;
  int cavitymode;
  int fbmode;
  int buffersize;
  int windowlength;
  double normfact;
  double tunergain;
  double *gain;
  double *turnhistory;
  double *z_cuts;
  double Length;
  double Voltage;
  double Energy;
  double Frequency;
  double HarmNumber;
  double TimeLag;
  double Qfactor;
  double Rshunt;
  double Beta;
  double phis;
  double feedback_angle_offset;
  double *vbunch;
  double *vbeam_phasor;
  double *vbeam;
  double *vcav;
  double *vgen;
  double *vgen_buffer;
  double *vbeam_buffer;
  double *vbunch_buffer;
  int system_harmonic;
  double ts;
  int every;
  int delay;
  int samplenum;
  int ff;
  double cutoff;
  int recordsize;
  double *Ig2Vg_vec;
  double *Ig2Vg_tmp;
  double *ig_phasor;
  double *ig_phasor_record;
  double *dot_output;
  double *generator_phasor_record;
  double *beam_phasor_record;
  double *cavity_phasor_record;
  double *Ig2Vg_mat;
  double *vc_previous;
  double *diff_record;
  double *samplelist;
  double *vc_list;
  double *I_record;
  double *FFconst;
  double *IIRout;
  double *IIRcoef;
  }; 


void write_buffer(double *data, double *buffer, int datasize, int buffersize){
    if(buffersize>1){
        memmove(buffer, buffer + datasize, datasize*(buffersize-1)*sizeof(double));
    }
    memcpy(buffer + datasize*(buffersize-1), data, datasize*sizeof(double));
}
   

void BeamLoadingCavityPass(double *r_in, int num_particles, int nbunch,
                           double *bunch_spos, double *bunch_currents, 
                           double *fillpattern,
                           double circumference,
                           int nturn, double energy, int harmonic_number,
                           int iturn,
                           struct elem *Elem) {
  
    long cavitymode = Elem->cavitymode;
    long fbmode = Elem->fbmode;
    
    long nslice = Elem->nslice;
    long nturnsw = Elem->nturnsw; /* can this attribute be removed? */
    long buffersize = Elem->buffersize;
    long windowlength = Elem->windowlength;

    double normfact = Elem->normfact;  
    double le = Elem->Length;
    double rffreq = Elem->Frequency;
    double harmn = Elem->HarmNumber;
    int ring_harmn = harmonic_number;
    double tlag = Elem->TimeLag;
    double qfactor = Elem->Qfactor;
    double rshunt = Elem->Rshunt;
    double beta = Elem->Beta;
    double tunergain = Elem->tunergain;
    //if fb mode is PROP then gain[0] is Voltgain and gain[1] is PhaseGain
    //if fb mode is PROP_INTEGRAL then gain[0] is Prop gain and gain[1] is integral gain
    double *gain = Elem->gain;
    double ts = Elem->ts;
    double *vgen_arr = Elem->vgen; /* [vgen, thetag, psi, vgr] */
        
    double *turnhistory = Elem->turnhistory;
    double *vgen_buffer = Elem->vgen_buffer;
    double *vbeam_buffer = Elem->vbeam_buffer;
    double *vbunch_buffer = Elem->vbunch_buffer;
    
    
    size_t sztmp1 = sizeof(double)*buffersize*3; //vgen, theta_g, psi
    void *set_params = atMalloc(sztmp1); // This is a buffer of the actual params to set. Good for PID
    
    double cutoff = Elem->cutoff;  
    int delay = Elem->delay; 
    int every = Elem->every; 
    int FF = Elem->ff; 
    int samplenum = Elem->samplenum; 
    int record_size = Elem->recordsize;
    int samplelist_length = ring_harmn/every;

    /* Here we have to declare empty pointers for the PI Loop
    They have to be defined outside of an if statement*/
    
    double *Ig2Vg_vec_real = Elem->Ig2Vg_vec;
    double *Ig2Vg_vec_imag = Elem->Ig2Vg_vec + ring_harmn;
    double *Ig2Vg_tmp_real = Elem->Ig2Vg_tmp;
    double *Ig2Vg_tmp_imag = Elem->Ig2Vg_tmp + ring_harmn;
    double *ig_phasor_real = Elem->ig_phasor;
    double *ig_phasor_imag = Elem->ig_phasor + ring_harmn;
    double *ig_phasor_record_real = Elem->ig_phasor_record;
    double *ig_phasor_record_imag = Elem->ig_phasor_record + ring_harmn;
    double *dot_output_real = Elem->dot_output;
    double *dot_output_imag = Elem->dot_output + ring_harmn;
                      

    double *generator_phasor_record_real = Elem->generator_phasor_record;
    double *generator_phasor_record_imag = Elem->generator_phasor_record + ring_harmn;
    double *beam_phasor_record_real = Elem->beam_phasor_record;
    double *beam_phasor_record_imag = Elem->beam_phasor_record + ring_harmn;
    double *cavity_phasor_record_real = Elem->cavity_phasor_record;
    double *cavity_phasor_record_imag = Elem->cavity_phasor_record + ring_harmn;

    double *Ig2Vg_mat_real = Elem->Ig2Vg_mat;
    double *Ig2Vg_mat_imag = Elem->Ig2Vg_mat + ring_harmn*ring_harmn;

    
    double *vc_previous_real = Elem->vc_previous; 
    double *vc_previous_imag = Elem->vc_previous + samplenum;
    double *diff_record_real = Elem->diff_record; 
    double *diff_record_imag = Elem->diff_record + record_size;
    double *samplelist = Elem->samplelist;
    double *vc_list_real = Elem->vc_list; 
    double *vc_list_imag = Elem->vc_list + ring_harmn + samplenum;
    
    //double *Ig_modulation_signal_real; double *Ig_modulation_signal_imag; 

    double *I_record = Elem->I_record;            
    double *FFconst = Elem->FFconst;
    double *IIRout = Elem->IIRout;
    double *IIRcoef = Elem->IIRcoef;

    double *z_cuts = Elem->z_cuts;
    double *vbunch = Elem->vbunch;
    double *vbeam = Elem->vbeam;
    double *vcav_set = Elem->vcav; /* Vcav set points amplitude, phase */
    double *vbeam_phasor = Elem->vbeam_phasor;
    double feedback_angle_offset = Elem->feedback_angle_offset;
    
    
        
        
    double vbeam_set[] = {vbeam[0], vbeam[1]};
    double vcav_meas[] = {0.0, 0.0, 0.0};
    double ave_vbeam[] = {0.0, 0.0};
    double tot_current = 0.0;

    
    int i;
    size_t sz = nslice*nbunch*sizeof(double) + num_particles*sizeof(int);
    int c;
    int *pslice;
    double *vbeam_kicks; /* This used to be kz, it is the kick that is applied */

    
    double vgen = vgen_arr[0];
    double gen_phase = vgen_arr[1];
    double psi = vgen_arr[2];
    double delta = pow(rffreq * tan(psi) / qfactor, 2) + 4 * pow(rffreq,2);
    double freqres = (rffreq * tan(psi) / qfactor + sqrt(delta)) / 2;

    double tot_lag_phase = (tlag+ts)*rffreq*TWOPI/C0;
    double filling_time = 2*qfactor / (TWOPI * freqres);
    double T1 = 1/rffreq;
    double kloss = rshunt * TWOPI * freqres / (2 * qfactor);

    double vcav_phasor[] = {0.0, 0.0}; 
    set_cavity_phasor(vgen, gen_phase, vbeam_phasor, vcav_phasor);
    
    for(i=0;i<nbunch;i++){
        tot_current += bunch_currents[i];
    }
    

    
    /*Track RF cavity is always done. */
    trackRFCavity(r_in, le, vgen/energy, rffreq, harmn, tlag, -gen_phase - tot_lag_phase, nturn, circumference/C0, num_particles);
    /*Only allocate memory if current is > 0*/
    if(tot_current>0){
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
                             beta, ave_vbeam, vbunch, bunch_spos, ring_harmn, fillpattern, ts);                        

        /*apply kicks*/
        for (c=0; c<num_particles; c++) {
            double *r6 = r_in+c*6;
            int islice=pslice[c];
            if (!atIsNaN(r6[0])) {         
                r6[4] += vbeam_kicks[islice];                 
            }
        }
        
               
        // First write the values to the buffer
        if(buffersize>0){
            write_buffer(vbeam, vbeam_buffer, 2, buffersize);
            write_buffer(vgen_arr, vgen_buffer, 4, buffersize);
            write_buffer(vbunch, vbunch_buffer, 2*ring_harmn, buffersize);
        }   

        update_vbeam_set(fbmode, vbeam_set, ave_vbeam, vbeam_buffer,
                             buffersize, windowlength);
                             
        compute_set_params(vbeam_set, vgen_arr, feedback_angle_offset, vcav_set[1], vcav_meas);
               
        if(cavitymode==1){
            // If CavityMode=ACTIVE
            if(tunergain>0){
                vgen_arr[2] += tunergain * (vcav_meas[2] - vgen_arr[2]);
            }
            if(fbmode==1){
                // If FBMode=PROP
                update_vgen(vcav_set, vgen_arr, vcav_meas, gain[0], gain[1], feedback_angle_offset);
            }
            if(fbmode==2){
                if(iturn==0){
                    init_sample_list(samplelist, ring_harmn, every); 
                             
                    init_phasor_arrays(vgen, gen_phase,
                                       ig_phasor_real, ig_phasor_imag,
                                       ig_phasor_record_real, ig_phasor_record_imag,
                                       ring_harmn, rshunt, psi,
                                       generator_phasor_record_real, generator_phasor_record_imag);

                    init_IIR(cutoff, IIRcoef, IIRout, T1, every, vcav_set[0]);
                    
                    init_FFconst(FF,
                                 ig_phasor_real, ig_phasor_imag,
                                 ring_harmn, FFconst);
                                 
                    init_Ig2Vg_matrix(ring_harmn,
                                      Ig2Vg_vec_real, Ig2Vg_vec_imag,
                                      Ig2Vg_tmp_real, Ig2Vg_tmp_imag,
                                      filling_time, psi, T1, 
                                      Ig2Vg_mat_real, Ig2Vg_mat_imag);

                    I_record[0] = 0.0; I_record[1] = 0.0;


                    set_cavity_phasor(vgen, gen_phase, ave_vbeam, vcav_phasor);

                    init_vc_previous(vc_previous_real, vc_previous_imag, samplenum, vcav_phasor);    

                };
                
                init_cavity_record_phasor_array(vbunch,
                                                beam_phasor_record_real, beam_phasor_record_imag,
                                                cavity_phasor_record_real, cavity_phasor_record_imag,
                                                generator_phasor_record_real, generator_phasor_record_imag,
                                                ring_harmn); 
                
                if(iturn>=1 && tunergain>0){
                    init_Ig2Vg_matrix(ring_harmn,
                                      Ig2Vg_vec_real, Ig2Vg_vec_imag,
                                      Ig2Vg_tmp_real, Ig2Vg_tmp_imag,
                                      filling_time, psi, T1,
                                      Ig2Vg_mat_real, Ig2Vg_mat_imag);
                }
                
                track_PIL(vc_previous_real, vc_previous_imag,
                          cavity_phasor_record_real, cavity_phasor_record_imag,
                          ig_phasor_real, ig_phasor_imag,
                          samplelist, samplenum, record_size, samplelist_length,
                          diff_record_real, diff_record_imag,
                          FFconst, gain, I_record,
                          rffreq,
                          vcav_set[0], vcav_set[1],
                          generator_phasor_record_real, generator_phasor_record_imag,
                          Ig2Vg_vec_real, Ig2Vg_vec_imag,
                          Ig2Vg_mat_real, Ig2Vg_mat_imag,
                          ig_phasor_record_real, ig_phasor_record_imag,
                          dot_output_real, dot_output_imag,
                          kloss, T1, ring_harmn, vgen_arr,
                          IIRout, IIRcoef,
                          vc_list_real, vc_list_imag,
                          every,
                          psi, rshunt
                          );    
                
            }
        }else if(cavitymode==3){     
            update_passive_frequency(vbeam_set, vcav_set, vgen_arr, tunergain);
        }

        vbeam[0] = ave_vbeam[0];
        vbeam[1] = ave_vbeam[1];
        atFree(buffer);
    }
    
}


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double rl = Param->RingLength;
    double energy;
    int nturn = Param->nturn;
    if (!Elem) {
        long nslice, nturns, cavitymode, fbmode, buffersize, windowlength, system_harmonic;
        long delay, every, samplenum, ff, recordsize;
        double wakefact, Energy, Frequency, TimeLag, Length, feedback_angle_offset;
        double normfact, tunergain, qfactor, rshunt, beta, phis, ts, cutoff;
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

        double *Ig2Vg_vec;
        double *Ig2Vg_tmp;
        double *ig_phasor;
        double *ig_phasor_record;
        double *dot_output;
        double *generator_phasor_record;
        double *beam_phasor_record;
        double *cavity_phasor_record;
        double *Ig2Vg_mat;
        double *vc_previous;
        double *diff_record;
        double *samplelist;
        double *vc_list;
        double *I_record;
        double *FFconst;
        double *IIRcoef;
        double *IIRout;
        
        /*attributes for RF cavity*/
        Length=atGetDouble(ElemData,"Length"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        /*attributes for resonator*/
        nslice=atGetLong(ElemData,"_nslice"); check_error();
        nturns=atGetLong(ElemData,"_nturns"); check_error();
        buffersize=atGetLong(ElemData,"_buffersize"); check_error();
        windowlength=atGetLong(ElemData,"_windowlength"); check_error();
        cavitymode=atGetLong(ElemData,"_cavitymode"); check_error();
        fbmode=atGetLong(ElemData,"_fbmode"); check_error();
        wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
        qfactor=atGetDouble(ElemData,"Qfactor"); check_error();
        rshunt=atGetDouble(ElemData,"Rshunt"); check_error();
        beta=atGetDouble(ElemData,"_beta"); check_error();
        normfact=atGetDouble(ElemData,"NormFact"); check_error();
        
        gain=atGetDoubleArray(ElemData,"Gain"); check_error();
        tunergain=atGetDouble(ElemData,"TunerGain"); check_error();
        
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
        system_harmonic=atGetLong(ElemData,"system_harmonic"); check_error();
        ts=atGetDouble(ElemData,"_ts"); check_error();
        
        /*optional attributes*/        
        delay=atGetOptionalLong(ElemData,"delay",1); check_error();
        every=atGetOptionalLong(ElemData,"every",1); check_error();
        samplenum=atGetOptionalLong(ElemData,"samplenum",1); check_error();        
        cutoff=atGetOptionalDouble(ElemData,"cutoff",0); check_error();        
        ff=atGetOptionalLong(ElemData,"FF",1); check_error();
        recordsize=atGetOptionalLong(ElemData,"recordsize",1); check_error();

        Ig2Vg_vec=atGetOptionalDoubleArray(ElemData,"_Ig2Vg_vec"); check_error();
        Ig2Vg_tmp=atGetOptionalDoubleArray(ElemData,"_Ig2Vg_tmp"); check_error();
        ig_phasor=atGetOptionalDoubleArray(ElemData,"_ig_phasor"); check_error();
        ig_phasor_record=atGetOptionalDoubleArray(ElemData,"_ig_phasor_record"); check_error();
        dot_output=atGetOptionalDoubleArray(ElemData,"_dot_output"); check_error();
        generator_phasor_record=atGetOptionalDoubleArray(ElemData,"_generator_phasor_record"); check_error();
        beam_phasor_record=atGetOptionalDoubleArray(ElemData,"_beam_phasor_record"); check_error();
        cavity_phasor_record=atGetOptionalDoubleArray(ElemData,"_cavity_phasor_record"); check_error();

        Ig2Vg_mat=atGetOptionalDoubleArray(ElemData,"_Ig2Vg_mat"); check_error();
        vc_previous=atGetOptionalDoubleArray(ElemData,"_vc_previous"); check_error();
        diff_record=atGetOptionalDoubleArray(ElemData,"_diff_record"); check_error();        
        samplelist=atGetOptionalDoubleArray(ElemData,"_samplelist"); check_error();        
        vc_list=atGetOptionalDoubleArray(ElemData,"_vc_list"); check_error();        
        I_record=atGetOptionalDoubleArray(ElemData,"_I_record"); check_error();
        FFconst=atGetOptionalDoubleArray(ElemData,"_FFconst"); check_error();
        IIRcoef=atGetOptionalDoubleArray(ElemData,"_IIRcoef"); check_error();
        IIRout=atGetOptionalDoubleArray(ElemData,"_IIRout"); check_error();
        
                
        Energy=atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
        feedback_angle_offset=atGetOptionalDouble(ElemData,"feedback_angle_offset", 0.0); check_error();

        /* Check energy */
        Energy = atEnergy(Param->energy, Energy); check_error();

        int dimsth[] = {Param->nbunch*nslice*nturns, 4};
        atCheckArrayDims(ElemData,"_turnhistory", 2, dimsth); check_error();
        int dimsvb[] = {Param->harmonic_number, 2};
        atCheckArrayDims(ElemData,"_vbunch", 2, dimsvb); check_error();
       
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        
        Elem->Length=Length;
        Elem->Frequency=Frequency;
        Elem->HarmNumber=round(Frequency*rl/C0);
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
        Elem->tunergain = tunergain;
        Elem->gain = gain;
        Elem->vbeam_phasor = vbeam_phasor;
        Elem->cavitymode = cavitymode;
        Elem->buffersize = buffersize;
        Elem->windowlength = windowlength;
        Elem->vgen_buffer = vgen_buffer;
        Elem->vbeam_buffer = vbeam_buffer;
        Elem->vbunch_buffer = vbunch_buffer;
        Elem->feedback_angle_offset = feedback_angle_offset;
        Elem->fbmode = fbmode;
        Elem->phis = phis;
        Elem->ts = ts;
        Elem->system_harmonic = system_harmonic;
        Elem->every=every;
        Elem->delay=delay;
        Elem->samplenum=samplenum;
        Elem->cutoff=cutoff;
        Elem->ff=ff;
        Elem->recordsize=recordsize;
        Elem->I_record=I_record;
        Elem->Ig2Vg_vec=Ig2Vg_vec;
        Elem->Ig2Vg_tmp=Ig2Vg_tmp;
        Elem->ig_phasor=ig_phasor;
        Elem->ig_phasor_record=ig_phasor_record;
        Elem->dot_output=dot_output;
        Elem->generator_phasor_record=generator_phasor_record;
        Elem->beam_phasor_record=beam_phasor_record;
        Elem->cavity_phasor_record=cavity_phasor_record;
        Elem->Ig2Vg_mat=Ig2Vg_mat;
        Elem->vc_previous=vc_previous;
        Elem->diff_record=diff_record;
        Elem->samplelist=samplelist;
        Elem->vc_list=vc_list;
        Elem->FFconst=FFconst;
        Elem->IIRcoef=IIRcoef;
        Elem->IIRout=IIRout;
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

      long nslice, nturns, cavitymode, fbmode, buffersize, windowlength, system_harmonic;
      long delay, every, samplenum, ff, recordsize;
      double wakefact, Energy, Frequency, TimeLag, Length, feedback_angle_offset;
      double normfact, tunergain, qfactor, rshunt, beta, phis, ts, cutoff;
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
      
      
      /*attributes for RF cavity*/
      Length=atGetDouble(ElemData,"Length"); check_error();
      Frequency=atGetDouble(ElemData,"Frequency"); check_error();
      TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
      /*attributes for resonator*/
      nslice=atGetLong(ElemData,"_nslice"); check_error();
      nturns=atGetLong(ElemData,"_nturns"); check_error();
      buffersize=atGetLong(ElemData,"_buffersize"); check_error();
      windowlength=atGetLong(ElemData,"_windowlength"); check_error();
      cavitymode=atGetLong(ElemData,"_cavitymode"); check_error();
      fbmode=atGetLong(ElemData,"_fbmode"); check_error();
      wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
      qfactor=atGetDouble(ElemData,"Qfactor"); check_error();
      rshunt=atGetDouble(ElemData,"Rshunt"); check_error();
      beta=atGetDouble(ElemData,"_beta"); check_error();
      normfact=atGetDouble(ElemData,"NormFact"); check_error();
      phasegain=atGetDouble(ElemData,"PhaseGain"); check_error();
      gain=atGetDouble(ElemData,"Gain"); check_error();
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
      system_harmonic=atGetLong(ElemData,"system_harmonic"); check_error();
      ts=atGetDouble(ElemData,"_ts"); check_error();

      /*optional attributes*/
      Energy=atGetOptionalDouble(ElemData,"Energy",0.0); check_error();
      z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
      feedback_angle_offset=atGetOptionalDouble(ElemData,"feedback_angle_offset",0.0); check_error();
      
      
      Ig2Vg_vec=atGetOptionalDoubleArray(ElemData,"_Ig2Vg_vec"); check_error();
      Ig2Vg_tmp=atGetOptionalDoubleArray(ElemData,"_Ig2Vg_tmp"); check_error();
      ig_phasor=atGetOptionalDoubleArray(ElemData,"_ig_phasor"); check_error();
      ig_phasor_record=atGetOptionalDoubleArray(ElemData,"_ig_phasor_record"); check_error();
      dot_output=atGetOptionalDoubleArray(ElemData,"_dot_output"); check_error();
      generator_phasor_record=atGetOptionalDoubleArray(ElemData,"_generator_phasor_record"); check_error();
      beam_phasor_record=atGetOptionalDoubleArray(ElemData,"_beam_phasor_record"); check_error();
      cavity_phasor_record=atGetOptionalDoubleArray(ElemData,"_cavity_phasor_record"); check_error();


      delay=atGetOptionalLong(ElemData,"delay",1); check_error();
      every=atGetOptionalLong(ElemData,"every",1); check_error();
      samplenum=atGetOptionalLong(ElemData,"sample_num",1); check_error();        
      cutoff=atGetOptionalDouble(ElemData,"cutoff",0); check_error();        
      ff=atGetOptionalLong(ElemData,"FF",1); check_error();
      recordsize=atGetOptionalLong(ElemData,"recordsize",1); check_error();  
        
      Ig2Vg_mat=atGetOptionalDoubleArray(ElemData,"_Ig2Vg_mat"); check_error();
      vc_previous=atGetOptionalDoubleArray(ElemData,"_vc_previous"); check_error();
      diff_record=atGetOptionalDoubleArray(ElemData,"_diff_record"); check_error();        
      samplelist=atGetOptionalIntArray(ElemData,"_samplelist"); check_error();        
      vc_list=atGetOptionalDoubleArray(ElemData,"_vc_list"); check_error();        
      I_record=atGetOptionalDoubleArray(ElemData,"_I_record"); check_error();
      FFconst=atGetOptionalDoubleArray(ElemData,"_FFcont"); check_error();
      IIRcoef=atGetOptionalDoubleArray(ElemData,"_IIRcoef"); check_error()
      IIRout=atGetOptionalDoubleArray(ElemData,"_IIRout"); check_error();
      
      Elem = (struct elem*)atMalloc(sizeof(struct elem));
      Elem->Length=Length;
      Elem->cavitymode=cavitymode;
      Elem->Frequency=Frequency;
      Elem->HarmNumber=1;
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
      Elem->tunergain = tunergain;
      Elem->vbeam_phasor = vbeam_phasor;
      Elem->buffersize = buffersize;
      Elem->windowlength = windowlength;
      Elem->vgen_buffer = vgen_buffer;
      Elem->vbeam_buffer = vbeam_buffer;
      Elem->vbunch_buffer = vbunch_buffer;
      Elem->feedback_angle_offset = feedback_angle_offset;
      Elem->phis = phis;
      Elem->ts = ts;
      Elem->system_harmonic = system_harmonic;
      Elem->every=every;
      Elem->delay=delay;
      Elem->samplenum=samplenum;
      Elem->cutoff=cutoff;
      Elem->ff=ff;
      Elem->recordsize=recordsize;
      
      Elem->Ig2Vg_vec=Ig2Vg_vec;
      Elem->Ig2Vg_tmp=Ig2Vg_tmp;
      Elem->ig_phasor=ig_phasor;
      Elem->ig_phasor_record=ig_phasor_record;
      Elem->dot_output=dot_output;
      Elem->generator_phasor_record=generator_phasor_record;
      Elem->beam_phasor_record=beam_phasor_record;
      Elem->cavity_phasor_record=cavity_phasor_record;
      Elem->Ig2Vg_mat=Ig2Vg_mat;
      Elem->vc_previous=vc_previous;
      Elem->diff_record=diff_record;
      Elem->samplelist=samplelist;
      Elem->vc_list=vc_list;
      Elem->I_record=I_record;
      Elem->FFconst=FFconst;
      Elem->IIRcoef=IIRcoef;
      Elem->IIRout=IIRout;
      
      Elem->fbmode = fbmode;
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
      mxSetCell(plhs[0],13,mxCreateString("VoltGain"));
      mxSetCell(plhs[0],14,mxCreateString("_turnhistory"));
      mxSetCell(plhs[0],15,mxCreateString("_vbunch"));
      mxSetCell(plhs[0],16,mxCreateString("_vbeam"));
      mxSetCell(plhs[0],17,mxCreateString("_vcav"));
      mxSetCell(plhs[0],18,mxCreateString("_vgen"));
      mxSetCell(plhs[0],19,mxCreateString("_vbeam_phasor"));
      mxSetCell(plhs[0],20,mxCreateString("_vgen_buffer"));
      mxSetCell(plhs[0],21,mxCreateString("_vbeam_buffer"));
      mxSetCell(plhs[0],22,mxCreateString("_vbunch_buffer"));
      mxSetCell(plhs[0],23,mxCreateString("_buffersize"));
      mxSetCell(plhs[0],24,mxCreateString("_windowlength"));
      mxSetCell(plhs[0],25,mxCreateString("_phis"));
      mxSetCell(plhs[0],26,mxCreateString("system_harmonic")); 
      mxSetCell(plhs[0],27,mxCreateString("_ts"));     
                                
      if(nlhs>1) /* optional fields */
      {
          plhs[1] = mxCreateCellMatrix(26,1);
          mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
          mxSetCell(plhs[1],1,mxCreateString("ZCuts"));
          mxSetCell(plhs[1],2,mxCreateString("feedback_angle_offset"));

          mxSetCell(plhs[1],3,mxCreateString("Delay"));           
          mxSetCell(plhs[1],4,mxCreateString("Every"));           
          mxSetCell(plhs[1],5,mxCreateString("SampleNum"));            
          mxSetCell(plhs[1],6,mxCreateString("Delay"));  
          mxSetCell(plhs[1],7,mxCreateString("FF"));   
          
          mxSetCell(plhs[1],8,mxCreateString("Ig2Vg_vec"));   
          mxSetCell(plhs[1],9,mxCreateString("Ig2Vg_tmp"));   
          mxSetCell(plhs[1],10,mxCreateString("ig_phasor"));   
          mxSetCell(plhs[1],11,mxCreateString("ig_phasor_record"));   
          mxSetCell(plhs[1],12,mxCreateString("dot_output"));   
          mxSetCell(plhs[1],13,mxCreateString("generator_phasor_record"));   
          mxSetCell(plhs[1],14,mxCreateString("beam_phasor_record"));   
          mxSetCell(plhs[1],15,mxCreateString("cavity_phasor_record"));   
          mxSetCell(plhs[1],16,mxCreateString("Ig2Vg_mat"));   
          mxSetCell(plhs[1],17,mxCreateString("vc_previous"));   
          mxSetCell(plhs[1],18,mxCreateString("diff_record"));   
          mxSetCell(plhs[1],19,mxCreateString("samplelist"));   
          mxSetCell(plhs[1],20,mxCreateString("vc_list"));   
          mxSetCell(plhs[1],21,mxCreateString("I_record"));   
          mxSetCell(plhs[1],22,mxCreateString("FFconst"));   
          mxSetCell(plhs[1],23,mxCreateString("IIRcoef"));   
          mxSetCell(plhs[1],24,mxCreateString("IIRout"));   
          mxSetCell(plhs[1],25,mxCreateString("RecordSize"));   
      }
  }
  else
  {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
  }
  
}
#endif /* MATLAB_MEX_FILE */

