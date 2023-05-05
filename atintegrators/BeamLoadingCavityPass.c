#include "atelem.c"
#include "atimplib.c"
#include "attrackfunc.c"

/*
 * BeamLoadingCavity pass method by Simon White.  
 *
 */
 
#define C0     2.99792458e8 

struct elem
{
  int nslice;
  int nturnsw;
  int mode;
  double normfact;
  double phasegain;
  double voltgain;
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
  double *vbunch;
  double *vbeam_phasor;
  double *vbeam;
  double *vcav;
  double *vgen;
}; 
   

void BeamLoadingCavityPass(double *r_in,int num_particles,int nbunch,
                           double *bunch_spos,double *bunch_currents,
                           double circumference,int nturn,struct elem *Elem) {
    /*
     * r_in - 6-by-N matrix of initial conditions reshaped into
     * 1-d array of 6*N elements
     */   
    long nslice = Elem->nslice;
    long nturnsw = Elem->nturnsw;
    long mode = Elem->mode;
    double normfact = Elem->normfact;  
    double le = Elem->Length;
    double energy = Elem->Energy;
    double rffreq = Elem->Frequency;
    double harmn = Elem->HarmNumber;
    double tlag = Elem->TimeLag;
    double qfactor = Elem->Qfactor;
    double rshunt = Elem->Rshunt;
    double beta = Elem->Beta;
    double *turnhistory = Elem->turnhistory;
    double *z_cuts = Elem->z_cuts;
    double *vbunch = Elem->vbunch;
    double phasegain = Elem->phasegain;
    double voltgain = Elem->voltgain;
    double *vbeam_phasor = Elem->vbeam_phasor;
    double *vbeamk = Elem->vbeam;
    double *vcavk = Elem->vcav;
    double *vgenk = Elem->vgen;
    double tot_current=0.0;
    int i;
    size_t sz = nslice*nbunch*sizeof(double) + num_particles*sizeof(int);
    int c;
    int *pslice;
    double *kz;
    double freqres = rffreq/(1-tan(vgenk[1])/(2*qfactor));
    double vgen = vgenk[0];
    double psi = vgenk[1];
    
    for(i=0;i<nbunch;i++){
        tot_current += bunch_currents[i];
    }

    /*Track RF cavity is always done. */
    trackRFCavity(r_in,le,vgen/energy,rffreq,harmn,tlag,-psi,nturn,circumference/C0,num_particles);

    
    /*Only allocate memory if current is > 0*/
    if(tot_current>0){
        void *buffer = atMalloc(sz);
        double *dptr = (double *) buffer;
        int *iptr;
        kz = dptr;
        dptr += nslice*nbunch;
        iptr = (int *) dptr;
        pslice = iptr; iptr += num_particles;
       
        rotate_table_history(nturnsw,nslice*nbunch,turnhistory,circumference);
        slice_bunch(r_in,num_particles,nslice,nturnsw,nbunch,bunch_spos,bunch_currents,
                    turnhistory,pslice,z_cuts);
        if(mode==2){
            compute_kicks_phasor(nslice,nbunch,nturnsw,turnhistory,normfact,kz,freqres,
                                 qfactor,rshunt,vbeam_phasor,circumference,energy,beta,
                                 vbeamk,vbunch);                        
        }else if(mode==1){
            compute_kicks_longres(nslice,nbunch,nturnsw,turnhistory,normfact,kz,freqres,
                                  qfactor,rshunt,beta,vbeamk,energy,vbunch);
        }
        /*apply kicks and RF*/
        /* OpenMP not efficient. Too much shared data ?
        #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
        shared(r_in,num_particles,pslice,kz) private(c)
        */   
        for (c=0; c<num_particles; c++) {
            double *r6 = r_in+c*6;
            int islice=pslice[c];
            if (!atIsNaN(r6[0])) {         
                r6[4] += kz[islice]; 
            }
        }
        update_vgen(vbeamk,vcavk,vgenk,phasegain,voltgain); 
        atFree(buffer);
    }
}


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double rl = Param->RingLength;
    int nturn=Param->nturn;
    if (!Elem) {
        long nslice,nturns,mode;
        double wakefact;
        double normfact, phasegain, voltgain;
        double *turnhistory;
        double *z_cuts;
        double Energy, Frequency, TimeLag, Length;
        double qfactor,rshunt,beta;
        double *vbunch;
        double *vbeam_phasor;
        double *vbeam;
        double *vgen;
        double *vcav;

        /*attributes for RF cavity*/
        Length=atGetDouble(ElemData,"Length"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        /*attributes for resonator*/
        nslice=atGetLong(ElemData,"_nslice"); check_error();
        nturns=atGetLong(ElemData,"_nturns"); check_error();
        mode=atGetLong(ElemData,"_mode"); check_error();
        wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
        qfactor=atGetDouble(ElemData,"Qfactor"); check_error();
        rshunt=atGetDouble(ElemData,"Rshunt"); check_error();
        beta=atGetDouble(ElemData,"_beta"); check_error();
        normfact=atGetDouble(ElemData,"NormFact"); check_error();
        phasegain=atGetDouble(ElemData,"PhaseGain"); check_error();
        voltgain=atGetDouble(ElemData,"VoltGain"); check_error();
        turnhistory=atGetDoubleArray(ElemData,"_turnhistory"); check_error();
        vbunch=atGetDoubleArray(ElemData,"_vbunch"); check_error();
        vbeam=atGetDoubleArray(ElemData,"_vbeam"); check_error();
        vcav=atGetDoubleArray(ElemData,"_vcav"); check_error();
        vgen=atGetDoubleArray(ElemData,"_vgen"); check_error();
        vbeam_phasor=atGetDoubleArray(ElemData,"_vbeam_phasor"); check_error(); 
        /*optional attributes*/
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
       
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        
        Elem->Length=Length;
        Elem->mode=mode;
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
        Elem->phasegain = phasegain;
        Elem->voltgain = voltgain;
        Elem->vbeam_phasor = vbeam_phasor;
    }
    if(num_particles<Param->nbunch){
        atError("Number of particles has to be greater or equal to the number of bunches.");
    }else if (num_particles%Param->nbunch!=0){
        atWarning("Number of particles not a multiple of the number of bunches: uneven bunch load.");
    }
    #ifdef _MSC_VER
    if(Elem->mode==2){
        atError("Beam loading Phasor mode not implemented in Windows.");
    }
    #endif
    BeamLoadingCavityPass(r_in,num_particles,Param->nbunch,Param->bunch_spos,
                          Param->bunch_currents,rl,nturn,Elem);
    return Elem;
}

MODULE_DEF(BeamLoadingCavityPass)       /* Dummy module initialisation */
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	
  if(nrhs == 2)
  {
  
      double *r_in;
      const mxArray *ElemData = prhs[0];
      int num_particles = mxGetN(prhs[1]);
      struct elem El, *Elem=&El;
      
      long nslice,nturns,mode;
      double wakefact;
      double normfact, phasegain, voltgain;
      double *turnhistory;
      double *z_cuts;
      double Energy, Frequency, TimeLag, Length;
      double qfactor,rshunt,beta;
      double *vbunch;
      double *vbeam_phasor;
      double *vbeam;
      double *vgen;
      double *vcav;
      /*attributes for RF cavity*/
      Length=atGetDouble(ElemData,"Length"); check_error();
      Energy=atGetDouble(ElemData,"Energy"); check_error();
      Frequency=atGetDouble(ElemData,"Frequency"); check_error();
      TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
      /*attributes for resonator*/
      nslice=atGetLong(ElemData,"_nslice"); check_error();
      nturns=atGetLong(ElemData,"_nturns"); check_error();
      mode=atGetLong(ElemData,"_mode"); check_error();
      wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
      qfactor=atGetDouble(ElemData,"Qfactor"); check_error();
      rshunt=atGetDouble(ElemData,"Rshunt"); check_error();
      beta=atGetDouble(ElemData,"_beta"); check_error();
      normfact=atGetDouble(ElemData,"NormFact"); check_error();
      phasegain=atGetDouble(ElemData,"PhaseGain"); check_error();
      voltgain=atGetDouble(ElemData,"VoltGain"); check_error();
      turnhistory=atGetDoubleArray(ElemData,"_turnhistory"); check_error();
      vbunch=atGetDoubleArray(ElemData,"_vbunch"); check_error();
      vbeam=atGetDoubleArray(ElemData,"_vbeam"); check_error();
      vcav=atGetDoubleArray(ElemData,"_vcav"); check_error();
      vgen=atGetDoubleArray(ElemData,"_vgen"); check_error();
      vbeam_phasor=atGetDoubleArray(ElemData,"_vbeam_phasor"); check_error(); 
      /*optional attributes*/
      z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
      Elem = (struct elem*)atMalloc(sizeof(struct elem));
      Elem->Length=Length;
      Elem->mode=mode;
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
      Elem->phasegain = phasegain;
      Elem->voltgain = voltgain;
      Elem->vbeam_phasor = vbeam_phasor;
      if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
      /* ALLOCATE memory for the output array of the same size as the input  */
      plhs[0] = mxDuplicateArray(prhs[1]);
      r_in = mxGetDoubles(plhs[0]);

      double bspos = 0.0;
      double bcurr = 0.0;
      BeamLoadingCavityPass(r_in,num_particles,1,&bspos,&bcurr,1,0,Elem);
  }
  else if (nrhs == 0)
  {   /* return list of required fields */
      plhs[0] = mxCreateCellMatrix(19,1);
      mxSetCell(plhs[0],0,mxCreateString("Length"));
      mxSetCell(plhs[0],1,mxCreateString("Energy"));
      mxSetCell(plhs[0],2,mxCreateString("Frequency"));
      mxSetCell(plhs[0],3,mxCreateString("_nslice"));
      mxSetCell(plhs[0],4,mxCreateString("_nturns"));
      mxSetCell(plhs[0],5,mxCreateString("_mode"));
      mxSetCell(plhs[0],6,mxCreateString("_wakefact"));
      mxSetCell(plhs[0],7,mxCreateString("Qfactor"));
      mxSetCell(plhs[0],8,mxCreateString("Rshunt"));
      mxSetCell(plhs[0],9,mxCreateString("_beta"));
      mxSetCell(plhs[0],10,mxCreateString("NormFact"));
      mxSetCell(plhs[0],11,mxCreateString("PhaseGain"));
      mxSetCell(plhs[0],12,mxCreateString("VoltGain"));
      mxSetCell(plhs[0],13,mxCreateString("_turnhistory"));
      mxSetCell(plhs[0],14,mxCreateString("_vbunch"));
      mxSetCell(plhs[0],15,mxCreateString("_vbeam"));
      mxSetCell(plhs[0],16,mxCreateString("_vcav"));
      mxSetCell(plhs[0],17,mxCreateString("_vgen"));
      mxSetCell(plhs[0],18,mxCreateString("_vbeam_phasor"));
      if(nlhs>1) /* optional fields */
      {
          plhs[1] = mxCreateCellMatrix(2,1);
          mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
          mxSetCell(plhs[0],1,mxCreateString("ZCuts"));
      }
  }
  else
  {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
  }
  
}
#endif /* MATLAB_MEX_FILE */
