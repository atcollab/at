/* AT interface to the elegant CCBEND deterministic orbital model, 2026.
 * See ccbend_core.h, elegant_ccbend_fringe.h and doc/licenseFiles/ElegantLicense.txt.
 * Adaptation: AT canonical coordinates, per-element optimization/cache,
 * normal/skew AT polynomials and AT alignment/aperture interfaces. */
#include "atelem.c"
#include "atlalib.c"
#include "ccbend_core.h"
#if CCB_RADIATION_MODE == 2
#include "atrandom.c"
typedef struct {pcg32_random_t *state; double cutoff;} CCBGaussian;
static double ccb_gaussian(void *ctx) {
    CCBGaussian *g=(CCBGaussian *)ctx;
    double sample;
    do {sample=atrandn_r(g->state,0,1);} while(fabs(sample)>g->cutoff);
    return sample;
}
#endif

#ifdef PYAT
#define ccb_check_error() if (PyErr_Occurred()) return 0
#else
#define ccb_check_error()
#endif

struct elem {
    CCBend bend;
    double Energy, GaussianCutoff;
    int QuantumSingleParticle;
    double *R1,*R2,*T1,*T2,*RApertures,*EApertures;
};

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
static double *ccb_array(const atElem *data, const char *name, int required,
                         int minimum) {
    int m=0,n=0;
    double *v=required?atGetDoubleArraySz(data,(char *)name,&m,&n):
        atGetOptionalDoubleArraySz(data,(char *)name,&m,&n);
    if(v && m*(n?n:1)<minimum) {
        atError("%s must contain at least %d doubles",name,minimum);
        return NULL;
    }
    return v;
}

static int ccb_read(const atElem *data, struct elem *out) {
    CCBend *e=&out->bend;
    memset(out,0,sizeof(*out));
    out->Energy=atGetOptionalDouble(data,"Energy",0);ccb_check_error();
    out->GaussianCutoff=atGetOptionalDouble(data,"GaussianCutoff",3);ccb_check_error();
    out->QuantumSingleParticle=atGetOptionalLong(data,"QuantumSingleParticle",1);ccb_check_error();
    if(!isfinite(out->GaussianCutoff) || out->GaussianCutoff<1) {
        atError("GaussianCutoff must be finite and at least 1");return 0;
    }
    e->length=atGetDouble(data,"Length");ccb_check_error();
    e->angle=atGetDouble(data,"BendingAngle");ccb_check_error();
    int maxorder=atGetLong(data,"MaxOrder");ccb_check_error();
    if(maxorder<0 || maxorder>8) {atError("MaxOrder must be between 0 and 8");return 0;}
    double *b=ccb_array(data,"PolynomB",1,maxorder+1);ccb_check_error();
    double *a=ccb_array(data,"PolynomA",1,maxorder+1);ccb_check_error();
    if(!a || !b) return 0;
    for(int i=0;i<=maxorder;i++) {e->b[i]=b[i];e->a[i]=a[i];}
    e->slices=atGetLong(data,"NumIntSteps");ccb_check_error();
    e->order=atGetOptionalLong(data,"IntegrationOrder",4);ccb_check_error();
    e->fringe=atGetOptionalLong(data,"FringeModel",0);ccb_check_error();
    e->edge_order=atGetOptionalLong(data,"EdgeOrder",3);ccb_check_error();
    e->edge_flip=atGetOptionalLong(data,"EdgeFlip",0);ccb_check_error();
    e->optimize_fse=atGetOptionalLong(data,"OptimizeFSE",1);ccb_check_error();
    e->optimize_dx=atGetOptionalLong(data,"OptimizeDX",1);ccb_check_error();
    e->compensate_kn=atGetOptionalLong(data,"CompensateKN",0);ccb_check_error();
    e->correction=atGetOptionalLong(data,"ReferenceCorrection",1);ccb_check_error();
    e->yaw=atGetOptionalDouble(data,"Yaw",0);ccb_check_error();
    e->gap=atGetOptionalDouble(data,"FullGap",0);ccb_check_error();
    e->fint[0]=atGetOptionalDouble(data,"FringeInt1",0);ccb_check_error();
    e->fint[1]=atGetOptionalDouble(data,"FringeInt2",0);ccb_check_error();
    e->fse=atGetOptionalDouble(data,"FSE",0);ccb_check_error();
    e->fse_dipole=atGetOptionalDouble(data,"FSEDipole",0);ccb_check_error();
    e->fse_quad=atGetOptionalDouble(data,"FSEQuadrupole",0);ccb_check_error();
    e->fse_offset=atGetOptionalDouble(data,"FSEOffset",0);ccb_check_error();
    e->dx_offset=atGetOptionalDouble(data,"DXOffset",0);ccb_check_error();
    /* Public order follows elegant parameter files: K0 K2 K4 K5 K6 K7 I0 I1.
     * Original elegant end-map routines use K0 I0 K2 I1 K4 K5 K6 K7. */
    static const int fringe_index[]={0,6,1,7,2,3,4,5};
    const char *fringe_names[]={"fringeInt1","fringeInt2"};
    for(int end=0;end<2;end++) {
        double *v=ccb_array(data,fringe_names[end],0,8);ccb_check_error();
        if(v) for(int i=0;i<8;i++) e->integrals[end][i]=v[fringe_index[i]];
    }
    out->R1=ccb_array(data,"R1",0,36);ccb_check_error();
    out->R2=ccb_array(data,"R2",0,36);ccb_check_error();
    out->T1=ccb_array(data,"T1",0,6);ccb_check_error();
    out->T2=ccb_array(data,"T2",0,6);ccb_check_error();
    out->RApertures=ccb_array(data,"RApertures",0,4);ccb_check_error();
    out->EApertures=ccb_array(data,"EApertures",0,2);ccb_check_error();
    if(!isfinite(e->length) || e->length<=0 || !isfinite(e->angle) || e->angle==0 ||
       fabs(e->angle)>=3.141592653589793 || e->slices<=0 ||
       (e->order!=2 && e->order!=4 && e->order!=6) || e->fringe < -1 || e->fringe>1 ||
       e->edge_order<0 || e->edge_order>3 || e->correction<0 || e->correction>3) {
        atError("Invalid CCBend: require Length>0, 0<abs(BendingAngle)<pi, NumIntSteps>0, "
                "IntegrationOrder=2/4/6, FringeModel=-1/0/1, EdgeOrder=0..3, ReferenceCorrection=0..3");
        return 0;
    }
    for(int i=0;i<9;i++) if(!isfinite(e->a[i]) || !isfinite(e->b[i])) {
        atError("PolynomA/B must be finite");return 0;
    }
    for(int j=0;j<2;j++) for(int i=0;i<8;i++) if(!isfinite(e->integrals[j][i])) {
        atError("fringeInt1/2 must be finite");return 0;
    }
    if(!ccb_prepare(e)) {atError("CCBend reference trajectory optimization failed");return 0;}
    return 1;
}

static int ccb_radiation(struct elem *e,double energy,double rest_energy,double charge) {
#if CCB_RADIATION_MODE > 0
    double gamma=atGamma(energy,e->Energy,rest_energy);ccb_check_error();
    if(!isfinite(gamma) || gamma<=1 || !isfinite(charge)) {
        atError("Radiating CCBEND requires finite Energy greater than rest energy");return 0;
    }
    double mass=rest_energy>0?rest_energy:__E0*1e9;
    double radius=__RE*(__E0*1e9/mass)*charge*charge;
    double po=sqrt(gamma*gamma-1);
    e->bend.po=po;
    e->bend.rad_coef=(2.0/3)*radius*po*po*po;
    /* This is the coefficient used by the local elegant ccbend.c ISR block.
     * AT CODATA constants are used for electron radius/rest energy. */
    e->bend.isr_coef=radius*sqrt(55.0/(24*sqrt(3.0))*pow(po,5)*137.0359895);
    if(e->bend.length<1e-6) e->bend.rad_coef=e->bend.isr_coef=0;
#endif
    return 1;
}
#endif

static void ccb_pass(double *r, int count, const struct elem *e) {
    for(int i=0;i<count;i++) {
        double *p=r+6*i;
        if(atIsNaN(p[0])) continue;
        if(e->T1) ATaddvv(p,e->T1);
        if(e->R1) ATmultmv(p,e->R1);
        if(e->RApertures) checkiflostRectangularAp(p,e->RApertures);
        if(e->EApertures) checkiflostEllipticalAp(p,e->EApertures);
        if(!atIsNaN(p[0]) && !ccb_at_track(p,&e->bend)) p[0]=atGetNaN();
        if(e->RApertures) checkiflostRectangularAp(p,e->RApertures);
        if(e->EApertures) checkiflostEllipticalAp(p,e->EApertures);
        if(e->R2) ATmultmv(p,e->R2);
        if(e->T2) ATaddvv(p,e->T2);
    }
}

typedef struct {int slice,next;double value;} CCBImpulse;
static double ccb_impulse(void *context) {
    CCBImpulse *imp=(CCBImpulse *)context;
    return imp->next++==imp->slice?imp->value:0;
}

/* Linear diffusion in AT coordinates. Source noise occurs at each slice end.
 * Differentiate the COMPLETE downstream map, including end maps, alignment
 * and radiation damping. This avoids replacing CCBEND by a sector-bend model.
 * Only used for envelope calculations; ordinary tracking has no extra cost. */
static int ccb_diffusion(const double *input,const struct elem *e,double *bdiff) {
    struct elem nominal=*e;
    nominal.bend.gaussian=NULL;nominal.bend.rng=NULL;
    double result[36]={0},jac[36],plus[6],minus[6];
    int incoming=0;
    for(int i=0;i<36;i++) if(bdiff[i]!=0) incoming=1;
    if(incoming) {
        for(int j=0;j<6;j++) {
            double h=j==5?1e-5:1e-7;
            memcpy(plus,input,sizeof(plus));memcpy(minus,input,sizeof(minus));
            plus[j]+=h;minus[j]-=h;
            ccb_pass(plus,1,&nominal);ccb_pass(minus,1,&nominal);
            if(!ccb_finite(plus) || !ccb_finite(minus)) return 0;
            for(int i=0;i<6;i++) jac[i+6*j]=(plus[i]-minus[i])/(2*h);
        }
        for(int i=0;i<6;i++) for(int j=0;j<6;j++)
            for(int k=0;k<6;k++) for(int l=0;l<6;l++)
                result[i+6*j]+=jac[i+6*k]*bdiff[k+6*l]*jac[j+6*l];
    }
    if(nominal.bend.isr_coef>0) {
        const double h=.01;
        CCBImpulse impulse={0,0,0};
        nominal.bend.gaussian=ccb_impulse;nominal.bend.rng=&impulse;
        for(int slice=0;slice<e->bend.slices;slice++) {
            memcpy(plus,input,sizeof(plus));memcpy(minus,input,sizeof(minus));
            impulse.slice=slice;impulse.next=0;impulse.value=h;
            ccb_pass(plus,1,&nominal);
            impulse.next=0;impulse.value=-h;ccb_pass(minus,1,&nominal);
            if(!ccb_finite(plus) || !ccb_finite(minus)) return 0;
            for(int i=0;i<6;i++) plus[i]=(plus[i]-minus[i])/(2*h);
            for(int i=0;i<6;i++) for(int j=0;j<6;j++)
                result[i+6*j]+=plus[i]*plus[j];
        }
    }
    memcpy(bdiff,result,sizeof(result));
    return 1;
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *data, struct elem *e,
        double *r, int count, struct parameters *param) {
    if(!e) {
        struct elem temp;
        if(!ccb_read(data,&temp)) return NULL;
        e=(struct elem *)atMalloc(sizeof(*e));*e=temp;
    }
    struct elem active=*e;
    if(!ccb_radiation(&active,param->energy,param->rest_energy,param->charge)) return NULL;
    if(param->bdiff) {
        if(count!=1 || !ccb_diffusion(r,&active,param->bdiff)) {
            atError("CCBEND diffusion requires one valid reference particle");return NULL;
        }
    }
#if CCB_RADIATION_MODE == 2
    CCBGaussian gaussian={param->thread_rng?param->thread_rng:&pcg32_global,active.GaussianCutoff};
    active.bend.gaussian=ccb_gaussian;active.bend.rng=&gaussian;
    if(count==1 && !active.QuantumSingleParticle) active.bend.gaussian=NULL;
#endif
    ccb_pass(r,count,&active);
    return e;
}
MODULE_DEF(CCB_MODULE)
#endif

#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    if(nrhs==2 || nrhs==3) {
        if(nlhs!=1 || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetM(prhs[1])!=6)
            mexErrMsgIdAndTxt("AT:WrongArg","Expected one output and a real double 6-by-N particle array");
        struct elem e;
        if(!ccb_read(prhs[0],&e)) return;
        double energy=0,rest_energy=0,charge=-1;
        if(nrhs==3) atProperties(prhs[2],&energy,&rest_energy,&charge);
        if(!ccb_radiation(&e,energy,rest_energy,charge)) return;
#if CCB_RADIATION_MODE == 2
        CCBGaussian gaussian={&pcg32_global,e.GaussianCutoff};
        e.bend.gaussian=ccb_gaussian;e.bend.rng=&gaussian;
        if(mxGetN(prhs[1])==1 && !e.QuantumSingleParticle) e.bend.gaussian=NULL;
#endif
        plhs[0]=mxDuplicateArray(prhs[1]);
        ccb_pass(mxGetDoubles(plhs[0]),(int)mxGetN(prhs[1]),&e);
    } else if(nrhs==0) {
        const char *req[]={"Length","BendingAngle","PolynomA","PolynomB","MaxOrder","NumIntSteps"};
        const char *opt[]={"IntegrationOrder","FringeModel","EdgeOrder","fringeInt1", "fringeInt2",
            "FringeInt1","FringeInt2","FullGap","EdgeFlip","Yaw","OptimizeFSE","OptimizeDX",
            "CompensateKN","ReferenceCorrection","FSE","FSEDipole","FSEQuadrupole","FSEOffset","DXOffset",
            "R1","R2","T1","T2","RApertures","EApertures", "Energy", "GaussianCutoff", "QuantumSingleParticle"};
        if(nlhs>0) {
            plhs[0]=mxCreateCellMatrix(6,1);
            for(int i=0;i<6;i++) mxSetCell(plhs[0],i,mxCreateString(req[i]));
        }
        if(nlhs>1) {
            int n=sizeof(opt)/sizeof(opt[0]);plhs[1]=mxCreateCellMatrix(n,1);
            for(int i=0;i<n;i++) mxSetCell(plhs[1],i,mxCreateString(opt[i]));
        }
    } else mexErrMsgIdAndTxt("AT:WrongArg","Needs 0, 2 or 3 arguments");
}
#endif
