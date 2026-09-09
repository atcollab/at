#include "atconstants.h"
#include "atelem.c"
#include "atimplib.c"
#include "attrackfunc.c"
#include "atfeedbacklib.c"
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
  double normfact;
  double phasegain;
  double voltgain;
  double *turnhistory;
  double *z_cuts;
  int delay;  double *VoltDelay; double *PhaseDelay;
  double *gain;
  double TunerOffset; int TunerAveragingPeriod; double TunerGain; double *TunerParams;  
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

  double *vbunch;
  double *vbeam_phasor;
  double *vbeam;
  double *vcav;
  double *vgen;
  double *vgen_buffer;
  double *vbeam_buffer;
  double *vbunch_buffer;
  double ts;
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

    double normfact = Elem->normfact;  
    double le = Elem->Length;
    double rffreq = Elem->Frequency;
    int harmn = rffreq * circumference / C0 ;    // cavity harmonic number 
       
    int ring_harmn = harmonic_number;
    double tlag = Elem->TimeLag;
    double qfactor = Elem->Qfactor;
    double rshunt = Elem->Rshunt;
    double beta = Elem->Beta; //not cavity beta
    
    //Tuner Variables
    double TunerGain = Elem->TunerGain;
    double TunerOffset = Elem->TunerOffset;
    int TunerAveragingPeriod = Elem->TunerAveragingPeriod;
    double *TunerParams = Elem->TunerParams; //TunerParams[0] is TunerCount, TunerParams[1] is TunerDiff

    //if fb mode is PROP then gain[0] is Voltgain and gain[1] is PhaseGain
    //if fb mode is PROP_INTEGRAL then gain[0] is Prop gain and gain[1] is integral gain
    double *gain = Elem->gain;

    int delay = Elem->delay; 
    double *VoltDelay = Elem->VoltDelay;
    double *PhaseDelay = Elem->PhaseDelay;

    double ts = Elem->ts;

    
    
    
    double *turnhistory = Elem->turnhistory;
    double *vgen_buffer = Elem->vgen_buffer;
    double *vbeam_buffer = Elem->vbeam_buffer;
    double *vbunch_buffer = Elem->vbunch_buffer;
    
    double *z_cuts = Elem->z_cuts;
    double *vbunch = Elem->vbunch;
    double *vbeam_phasor = Elem->vbeam_phasor;
    double *vbeam = Elem->vbeam;
    double *vcav_set = Elem->vcav; /* Vcav set points amplitude, phase */


    
    double vbeam_set[] = {vbeam[0], vbeam[1]};
    double vcav_meas[] = {0.0, 0.0, 0.0};
    double ave_vbeam[] = {0.0, 0.0};
    double tot_current = 0.0;
    
    int i;
    size_t sz = nslice*nbunch*sizeof(double) + num_particles*sizeof(int);
    int c;
    int *pslice;
    double *vbeam_kicks; /* This used to be kz, it is the kick that is applied */
    double *vgen_arr = Elem->vgen; /* [vgen, thetag, psi, vgr] */
    
    double vgen = vgen_arr[0];
    double gen_phase = vgen_arr[1];
    double psi = vgen_arr[2];
    double delta = pow(rffreq * tan(psi) / qfactor, 2) + 4 * pow(rffreq,2);
    double freqres = (rffreq * tan(psi) / qfactor + sqrt(delta)) / 2;

    double tot_lag_phase = (tlag+ts)*rffreq*TWOPI/C0;
    double filling_time = 2*qfactor / (TWOPI * freqres);
    double T1 = 1/rffreq;
    double kloss = rshunt * TWOPI * freqres / (2 * qfactor);

    for(i=0;i<nbunch;i++){
        tot_current += bunch_currents[i];
    }
        
        
    /*Track RF cavity is always done. */
    trackRFCavity(r_in, le, vgen/energy, rffreq, harmn, tlag, -gen_phase - tot_lag_phase, nturn, circumference/C0, num_particles);
    /*Only allocate memory if current is > 0*/
    if(tot_current>0 && rshunt>0){
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
            write_buffer(vgen_arr, vgen_buffer, 3, buffersize);
            write_buffer(vbunch, vbunch_buffer, 2*ring_harmn, buffersize);
        }   


        vbeam_set[0] = ave_vbeam[0];
        vbeam_set[1] = ave_vbeam[1];        

        compute_set_params(vbeam_set, vgen_arr, vcav_set[1], vcav_meas);
        
        if(cavitymode==1){
            update_vgen(vcav_set, vgen_arr, vcav_meas, gain[0], gain[1], VoltDelay, PhaseDelay, delay);

        }else if(cavitymode==3){     
            update_passive_frequency(vbeam_set, vcav_set, vgen_arr, TunerGain);
        }


        /* Here is where the tuner is calculated and applied */
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
    int nturn=Param->nturn;
    if (!Elem) {
        long nslice, nturns, cavitymode, fbmode, buffersize;
        double wakefact;
        double normfact;
        int delay;
        double TunerGain, TunerOffset, TunerAveragingPeriod, *TunerParams;
        double *VoltDelay, *PhaseDelay;
        double *gain;
        double *turnhistory;
        double *vgen_buffer;
        double *vbeam_buffer;
        double *vbunch_buffer;
        double *z_cuts;
        double Energy, Frequency, TimeLag, Length;
        double qfactor,rshunt,beta;
        double *vbunch;
        double *vbeam_phasor;
        double *vbeam;
        double *vgen;
        double *vcav;
        double phis;
        double ts;

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
        
        /*optional attributes*/
        delay=atGetOptionalLong(ElemData,"delay",1); check_error();        
        VoltDelay=atGetOptionalDoubleArray(ElemData,"VoltDelay"); check_error();
        PhaseDelay=atGetOptionalDoubleArray(ElemData,"PhaseDelay"); check_error();
        Energy=atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
        TunerOffset=atGetOptionalDouble(ElemData,"TunerOffset", 0.0); check_error();
        TunerAveragingPeriod=atGetOptionalDouble(ElemData,"TunerAveragingPeriod",1); check_error();
        TunerParams=atGetOptionalDoubleArray(ElemData,"_TunerParams"); check_error();

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
        Elem->z_cuts=z_cuts;
        Elem->vbunch = vbunch;
        Elem->vbeam = vbeam;
        Elem->vgen = vgen;
        Elem->vcav = vcav;
        Elem->TunerGain = TunerGain;
        Elem->TunerOffset = TunerOffset;
        Elem->TunerAveragingPeriod = TunerAveragingPeriod;
        Elem->TunerParams = TunerParams;
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
      
      long nslice, nturns, cavitymode, fbmode, buffersize, windowlength;
      long delay, every, samplenum, ff, recordsize, openloop;
      double TunerGain, TunerOffset, TunerAveragingPeriod, *TunerParams;
      double *VoltDelay, *PhaseDelay;
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
      
      /*optional attributes*/
      Energy=atGetOptionalDouble(ElemData,"Energy",0.0); check_error();
      z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
      TunerAveragingPeriod=atGetOptionalLong(ElemData,"TunerAveragingPeriod",1); check_error();
      TunerOffset=atGetOptionalDouble(ElemData,"TunerOffset",0.0); check_error();
      TunerParams=atGetOptionalDoubleArray(ElemData,"_TunerParams"); check_error();
      
      delay=atGetOptionalLong(ElemData,"delay",1); check_error();
      VoltDelay=atGetOptionalDoubleArray(ElemData,"VoltDelay"); check_error();
      PhaseDelay=atGetOptionalDoubleArray(ElemData,"PhaseDelay"); check_error();  
      
      
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
      plhs[0] = mxCreateCellMatrix(25,1);
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
      if(nlhs>1) /* optional fields */
      {
          plhs[1] = mxCreateCellMatrix(8,1);
          mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
          mxSetCell(plhs[1],1,mxCreateString("ZCuts"));
          mxSetCell(plhs[1],2,mxCreateString("TunerOffset"));
          mxSetCell(plhs[1],3,mxCreateString("TunerAveragingPeriod"));
          mxSetCell(plhs[1],4,mxCreateString("Delay"));     
          mxSetCell(plhs[1],5,mxCreateString("TunerParams"));   
          mxSetCell(plhs[1],6,mxCreateString("VoltDelay"));
          mxSetCell(plhs[1],7,mxCreateString("PhaseDelay"));           
      }
  }
  else
  {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
  }
  
}
#endif /* MATLAB_MEX_FILE */

