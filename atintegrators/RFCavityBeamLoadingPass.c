#include "atelem.c"
#include "atimplib.c"
#include "driftkickrad.c"
#include <math.h>
#include <float.h>
/*
 * RFCavityBeamLoading pass method by Simon White.  
 * User may contact simon.white@esrf.fr for questions and comments.
 */
 
#define TWOPI  4*asin(1)
#define C0     2.99792458e8 
#define QE     1.602176634e-19

struct elem
{
  int nslice;
  int nturnsw;
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
  double Phil;
  double Rshunt;
  double Beta;
  double phis;
  /*pointer [freqres, vbeam, vcav, vgen, psi]*/
  double *bl_params;
  double *vbunch;
};


void trackCavity(double *r_in, double le, double nv, double freq, double h, double lag, double philag,
                 int nturn, double T0, int num_particles) {
    int c;
    if (le == 0) {
        for (c = 0; c<num_particles; c++) {
            double *r6 = r_in+c*6;
            if(!atIsNaN(r6[0]))
                r6[4] += -nv*sin(TWOPI*freq*((r6[5]-lag)/C0 - (h/freq-T0)*nturn) - philag);
        }
    }
    else {
        double halflength = le/2;
        for (c = 0;c<num_particles;c++) {
            double *r6 = r_in+c*6;
            if(!atIsNaN(r6[0]))  {
                drift6(r6, halflength);
                r6[4] += -nv*sin(TWOPI*freq*((r6[5]-lag)/C0 - (h/freq-T0)*nturn) - philag);
                drift6(r6, halflength);
            }
        }
    }
}


double get_vbeam(int nslice, int nbunch, double energy, double *kz, double *vbunch){
    int i,ib;
    double vbeam=0.0;
    for(i=0; i<nbunch;i++){
        vbunch[i] = 0.0;
    }    
    for(i=0; i<nslice*nbunch;i++){
        ib = (int)(i/nslice);
        vbeam += kz[i]*energy;
        vbunch[ib] += kz[i]*energy;
    }
    return vbeam/nbunch;
}


void update_params(double vbeam, double qfactor, double rfv, double phis,
                   double phil, double rffreq, double *bl_params,
                   double *vbunch, int nbunch,double voltgain,double phasegain){ 
                                             
    double theta,a,b,vbn;
    double vgen,freqres,vcav,psi; 
    int i;
    psi = bl_params[4];
    for(i=0;i<nbunch;i++){
        vbunch[i] = -vbunch[i]/(cos(psi)*cos(psi));
    } 
    vbn  = -vbeam/(cos(psi)*cos(psi)); 
    theta = phis+phil;
    a = rfv*cos(phil);
    b = rfv*sin(phil)-vbn*cos(theta);
    psi = asin(b/sqrt(a*a+b*b));
    vgen = rfv*cos(psi)-vbeam/cos(psi)*sin(phis);
    bl_params[3] += (vgen-bl_params[3])*voltgain;
    bl_params[4] += (psi-bl_params[4])*phasegain;
    bl_params[0] = rffreq/(1-tan(bl_params[4])/(2*qfactor));
    bl_params[2] = (bl_params[3]*sin(theta-bl_params[4])+vbeam)/sin(phis);   
    bl_params[1] = vbn;
} 
   

void RFCavityBeamLoadingPass(double *r_in,int num_particles,int nbunch,
                             double *bunch_spos,double *bunch_currents,
                             double circumference,int nturn,struct elem *Elem) {
    /*
     * r_in - 6-by-N matrix of initial conditions reshaped into
     * 1-d array of 6*N elements
     */   
    long nslice = Elem->nslice;
    long nturnsw = Elem->nturnsw;
    double normfact = Elem->normfact;  
    double le = Elem->Length;
    double rfv = Elem->Voltage;
    double energy = Elem->Energy;
    double rffreq = Elem->Frequency;
    double harmn = Elem->HarmNumber;
    double tlag = Elem->TimeLag;
    double qfactor = Elem->Qfactor;
    double rshunt = Elem->Rshunt;
    double beta = Elem->Beta;
    double phis = Elem->phis;
    double phil = Elem->Phil;
    double *turnhistory = Elem->turnhistory;
    double *z_cuts = Elem->z_cuts;
    double *bl_params = Elem->bl_params;
    double *vbunch = Elem->vbunch;
    double phasegain = Elem->phasegain;
    double voltgain = Elem->voltgain;

    size_t sz = nslice*nbunch*sizeof(double) + num_particles*sizeof(int);
    int c;

    int *pslice;
    double *kz;
    double freqres = bl_params[0];
    double vbeam = bl_params[1];
    double vgen = bl_params[3];
    double psi = bl_params[4];
    
    trackCavity(r_in,le,vgen/energy,rffreq,harmn,tlag,-psi,nturn,circumference/C0,num_particles);
       
    void *buffer = atMalloc(sz);
    double *dptr = (double *) buffer;
    int *iptr;

    kz = dptr;
    dptr += nslice*nbunch;
    iptr = (int *) dptr;
    pslice = iptr; iptr += num_particles;

    /*slices beam and compute kick and update cavity params*/
    rotate_table_history(nturnsw,nslice*nbunch,turnhistory,circumference);
    slice_bunch(r_in,num_particles,nslice,nturnsw,nbunch,bunch_spos,bunch_currents,
                turnhistory,pslice,z_cuts);
    if(normfact!=0){
        compute_kicks_longres(nslice*nbunch,nturnsw,turnhistory,normfact,kz,freqres,qfactor,rshunt,beta);
        vbeam = get_vbeam(nslice,nbunch,energy,kz,vbunch);
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
    } else{
        vbeam = 0.0;
    }
    update_params(vbeam,qfactor,rfv,phis,phil,rffreq,bl_params,vbunch,nbunch,phasegain,voltgain);  
    atFree(buffer);
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
        double Voltage, Energy, Frequency, TimeLag, Length;
        double qfactor,rshunt,beta,phis,phil;
        double *bl_params;
        double *vbunch;

        /*attributes for RF cavity*/
        Length=atGetDouble(ElemData,"Length"); check_error();
        Voltage=atGetDouble(ElemData,"Voltage"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        /*attributes for resonator*/
        nslice=atGetLong(ElemData,"_nslice"); check_error();
        nturns=atGetLong(ElemData,"_nturns"); check_error();
        wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
        qfactor=atGetDouble(ElemData,"Qfactor"); check_error();
        rshunt=atGetDouble(ElemData,"Rshunt"); check_error();
        beta=atGetDouble(ElemData,"_beta"); check_error();
        phis=atGetDouble(ElemData,"_phis"); check_error();
        phil=atGetDouble(ElemData,"Phil"); check_error();
        normfact=atGetDouble(ElemData,"NormFact"); check_error();
        phasegain=atGetDouble(ElemData,"PhaseGain"); check_error();
        voltgain=atGetDouble(ElemData,"VoltGain"); check_error();
        turnhistory=atGetDoubleArray(ElemData,"_turnhistory"); check_error();
        bl_params=atGetDoubleArray(ElemData,"_bl_params"); check_error();
        vbunch=atGetDoubleArray(ElemData,"_vbunch"); check_error();
        /*optional attributes*/
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
       
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        
        Elem->Length=Length;
        Elem->Voltage=Voltage;
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
        Elem->Phil = phil;
        Elem->Beta = beta;
        Elem->phis = phis;
        Elem->bl_params = bl_params;
        Elem->z_cuts=z_cuts;
        Elem->vbunch = vbunch;
        Elem->phasegain = phasegain;
        Elem->voltgain = voltgain;
    }
    RFCavityBeamLoadingPass(r_in,num_particles,Param->nbunch,Param->bunch_spos,
                            Param->bunch_currents,rl,nturn,Elem);
    return Elem;
}

MODULE_DEF(RFCavityBeamLoadingPass)       /* Dummy module initialisation */
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/
