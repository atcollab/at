#include "atelem.c"
#include "atimplib.c"
#include "driftkickrad.c"
#include <math.h>
#include <float.h>
/*
 * RFCavityBeamLoading pass method by Simon White.  
 * User may contact simon.white@esrf.fr for questions and comments.
 */
 
#define TWOPI  6.28318530717959
#define C0     2.99792458e8 
#define QE     1.602176634e-19

struct elem
{
  int nslice;
  int nelem;
  int nturnsw;
  double normfact;
  double *turnhistory;
  double *z_cuts;
  double Length;
  double Voltage;
  double Energy;
  double Frequency;
  double HarmNumber;
  double TimeLag;
  double PhaseLag;
  double Qfactor;
  double Phil;
  double Rshunt;
  double Beta;
  double u0;
  double num_charges;
  /*pointer [freqres, vbeam, vcav, vgen, psi]*/
  double *bl_params;
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

double get_vbeam(int nslice, double normfact, double num_charge, double *kz, double psi){
    int i;
    double vbeam;
    vbeam=0;
    for(i=0; i<nslice;i++){
        vbeam += kz[i]/normfact/nslice;
    }
    vbeam = fabs(vbeam)*num_charge*QE/(cos(psi)*cos(psi));
    return vbeam;
}


void update_params(double vbeam, double qfactor, double rfv, double u0,
                   double phil, double rffreq, double *bl_params){ 
                                             
    double phis,theta,a,b,x,dff;
    double vgen,freqres,vcav,psi; 
     
    phis = asin(u0/rfv);
    theta = phis+phil;
    a = rfv*cos(phis)*cos(theta)+u0*sin(theta);
    b = rfv*cos(phis)*sin(theta)-(vbeam+u0)*cos(theta);
    x = a*b/sqrt(a*a*(a*a+b*b));
    psi = asin(x);
    vgen = rfv*cos(psi)+vbeam*cos(psi)*sin(phis);
    dff = -tan(psi)/(2*qfactor);
    freqres = rffreq/(dff+1);
    vcav = (vgen*sin(theta-psi)-vbeam*cos(psi)*cos(psi))/sin(theta);
    
    bl_params[0] = freqres;
    bl_params[1] = vbeam;
    bl_params[2] = vcav;
    bl_params[3] = vgen;
    bl_params[4] = psi;
}
     

void RFCavityBeamLoadingPass(double *r_in,int num_particles,double circumference,int nturn,struct elem *Elem) {
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
    double plag = Elem->PhaseLag;
    double qfactor = Elem->Qfactor;
    double rshunt = Elem->Rshunt;
    double beta = Elem->Beta;
    double u0 = Elem->u0;
    double phil = Elem->Phil;
    double num_charges = Elem->num_charges;
    double *turnhistory = Elem->turnhistory;
    double *z_cuts = Elem->z_cuts;
    double *bl_params = Elem->bl_params;

    size_t sz = nslice*sizeof(double) + num_particles*sizeof(int);
    int c;

    int *pslice;
    double *kz;
    double freqres = bl_params[0];
    double vbeam = bl_params[1];
    double vgen = bl_params[3];
    double psi = bl_params[4];
    
    if (plag) {
        trackCavity(r_in,le,rfv/energy,rffreq,harmn,tlag,-plag,nturn,circumference/C0,num_particles);
    } else {
        trackCavity(r_in,le,vgen/energy,rffreq,harmn,tlag,-psi,nturn,circumference/C0,num_particles);
    }

    void *buffer = atMalloc(sz);
    double *dptr = (double *) buffer;
    int *iptr;

    kz = dptr;
    dptr += nslice;
    iptr = (int *) dptr;
    pslice = iptr; iptr += num_particles;

    /*slices beam and compute kick and update cavity params*/
    rotate_table_history(nturnsw,nslice,turnhistory,circumference);
    slice_bunch(r_in,num_particles,nslice,nturnsw,turnhistory,pslice,z_cuts);
    compute_kicks_longres(nslice,nturnsw,turnhistory,normfact,kz,freqres,qfactor,rshunt,beta);
    vbeam = get_vbeam(nslice,normfact,num_charges,kz,psi);
    update_params(vbeam,qfactor,rfv,u0,phil,rffreq,bl_params);
   
    /*apply kicks and RF*/
    /* OpenMP not efficient. Too much shared data ?
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
    shared(r_in,num_particles,pslice,kx,kx2,ky,ky2,kz) private(c)
    */   
    for (c=0; c<num_particles; c++) {
        double *r6 = r_in+c*6;
        int islice=pslice[c];
        if (!atIsNaN(r6[0])) {         
            r6[4] += kz[islice]; 
        }
    }
    atFree(buffer);
}


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double rl = Param->RingLength;
    int nturn=Param->nturn;
    if (!Elem) {
        long nslice,nturns;
        double num_charges, wakefact;
        double normfact;
        double *turnhistory;
        double *z_cuts;
        double Voltage, Energy, Frequency, TimeLag, Length, PhaseLag;
        double qfactor,rshunt,beta,u0,phil;
        double *bl_params;

        /*attributes for RF cavity*/
        Length=atGetDouble(ElemData,"Length"); check_error();
        Voltage=atGetDouble(ElemData,"Voltage"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        PhaseLag=atGetOptionalDouble(ElemData,"PhaseLag",0); check_error();
        /*attributes for resonator*/
        nslice=atGetLong(ElemData,"_nslice"); check_error();
        nturns=atGetLong(ElemData,"_nturns"); check_error();
        num_charges=atGetDouble(ElemData,"NumParticles"); check_error();
        wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
        qfactor=atGetDouble(ElemData,"Qfactor"); check_error();
        rshunt=atGetDouble(ElemData,"Rshunt"); check_error();
        beta=atGetDouble(ElemData,"_beta"); check_error();
        u0=atGetDouble(ElemData,"_u0"); check_error();
        phil=atGetDouble(ElemData,"Phil"); check_error();
        normfact=atGetDouble(ElemData,"NormFact"); check_error();
        turnhistory=atGetDoubleArray(ElemData,"_turnhistory"); check_error();
        bl_params=atGetDoubleArray(ElemData,"_bl_params"); check_error();
        /*optional attributes*/
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
       
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        
        Elem->Length=Length;
        Elem->Voltage=Voltage;
        Elem->Frequency=Frequency;
        Elem->HarmNumber=round(Frequency*rl/C0);
        Elem->Energy = Energy;
        Elem->TimeLag=TimeLag;
        Elem->PhaseLag=PhaseLag;     
        Elem->nslice=nslice;
        Elem->nturnsw=nturns;
        Elem->normfact=normfact*num_charges*wakefact;
        Elem->num_charges = num_charges*normfact;
        Elem->turnhistory=turnhistory;
        Elem->Qfactor = qfactor;
        Elem->Rshunt = rshunt;
        Elem->Phil = phil;
        Elem->Beta = beta;
        Elem->u0 = u0;
        Elem->bl_params = bl_params;
        Elem->z_cuts=z_cuts;
    }
    RFCavityBeamLoadingPass(r_in,num_particles,rl,nturn,Elem);
    return Elem;
}

MODULE_DEF(RFCavityBeamLoadingPass)       /* Dummy module initialisation */
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/
