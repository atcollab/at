/* CCBEND adapter, 2026. See elegant_ccbend_fringe.h and elegant.LICENSE.
 * Portions adapted from elegant ccbend.c:
 * Copyright (c) 2018 The University of Chicago, as Operator of Argonne
 * National Laboratory, and The Regents of the University of California,
 * as Operator of Los Alamos National Laboratory.
 * Modified: AT polynomials, local reference optimization, AT coordinates,
 * callback-based Gaussian excitation, and independently assembled body map.
 * Internal coordinates are elegant [x,x',y,y',s,delta]. No global state.
 * Normal field coefficients are AT coefficients (Kn/n!). */
#ifndef AT_CCBEND_CORE_H
#define AT_CCBEND_CORE_H
#include <math.h>
#include <float.h>
#include <string.h>
#include "elegant_ccbend_fringe.h"

typedef struct {
    double length, angle, b[9], a[9], fint[2], gap, integrals[2][8];
    double yaw, fse, fse_dipole, fse_quad, fse_offset, dx_offset;
    double x_adjust, reference[5];
    double po, rad_coef, isr_coef;
    double (*gaussian)(void *);
    void *rng;
    int slices, order, fringe, edge_order, edge_flip;
    int optimize_fse, optimize_dx, compensate_kn, correction;
} CCBend;

static int ccb_finite(const double *r) {
    for (int i=0;i<6;i++) if (!isfinite(r[i])) return 0;
    return 1;
}

static void ccb_coefficients(const CCBend *e, double *b, double *a, int ideal) {
    double f=e->fse_offset+(ideal?0:e->fse);
    double factor=(1+f)/(e->compensate_kn?1+e->fse_offset:1);
    for (int i=0;i<9;i++) {
        b[i]=e->b[i]*factor; a[i]=e->a[i]*factor;
        if (e->angle<0) {
            if (!(i%2)) b[i]=-b[i];
            if (i%2) a[i]=-a[i];
        }
    }
    b[0]+=(1+f+(ideal?0:e->fse_dipole))*fabs(e->angle)/e->length;
    b[1]+=e->b[1]*(ideal?0:e->fse_quad)/(e->compensate_kn?1+e->fse_offset:1);
}

static int ccb_kick(double *r, double *px, double *py, const double *b,
                    const double *a, double ds) {
    double re=b[8], im=a[8];
    for(int i=7;i>=0;i--) {
        double old=re;
        re=old*r[0]-im*r[2]+b[i];
        im=im*r[0]+old*r[2]+a[i];
    }
    *px-=ds*re; *py+=ds*im;
    double z2=(1+r[5])*(1+r[5])-(*px)*(*px)-(*py)*(*py);
    if (!(z2>0) || !(1+r[5]>0)) return 0;
    r[1]=*px/sqrt(z2); r[3]=*py/sqrt(z2);
    return 1;
}

/* Mirrors integrate_kick_KnL: extrema are sampled at slice ends, as in elegant. */
static int ccb_body(double *r, const CCBend *e, const double *b, const double *a,
                    double *xmin, double *xmax, double *xfinal) {
    static const double d2[]={.5,.5}, k2[]={1,0};
    static const double d4[]={.6756035959798289,-.1756035959798288,
                              -.1756035959798288,.6756035959798289};
    static const double k4[]={1.3512071919596578,-1.7024143839193153,
                              1.3512071919596578,0};
    static const double d6[]={.39225680523878,.5100434119184585,
        -.47105338540975655,.0687531682525181,.0687531682525181,
        -.47105338540975655,.5100434119184585,.39225680523878};
    static const double k6[]={.784513610477560,.235573213359357,
        -1.17767998417887,1.3151863206839063,-1.17767998417887,
        .235573213359357,.784513610477560,0};
    const double *d=e->order==2?d2:(e->order==4?d4:d6);
    const double *k=e->order==2?k2:(e->order==4?k4:k6);
    int n=e->order==2?2:(e->order==4?4:8);
    double chord=2*e->length/e->angle*sin(e->angle/2), ds=chord/e->slices;
    double norm=(1+r[5])/sqrt(1+r[1]*r[1]+r[3]*r[3]);
    double px=norm*r[1], py=norm*r[3];
    double initial_s=r[4], path=0;
    double initial_p=e->po*(1+r[5]);
    double beta0=e->rad_coef?initial_p/sqrt(initial_p*initial_p+1):1;
    *xmin=DBL_MAX; *xmax=-DBL_MAX;
    for(int j=0;j<e->slices;j++) {
        double delta_px=0,delta_py=0;
        for(int i=0;i<n;i++) {
            double dl=ds*d[i];
            path+=dl*sqrt(1+r[1]*r[1]+r[3]*r[3]);
            r[0]+=dl*r[1]; r[2]+=dl*r[3];
            if (i<n-1) {
                double before_px=px,before_py=py;
                if(!ccb_kick(r,&px,&py,b,a,ds*k[i])) return 0;
                delta_px=px-before_px;delta_py=py-before_py;
            }
        }
        /* elegant integrate_kick_KnL applies radiation once per positive-length
         * slice, using the final kick's local field, not negative substeps. */
        if(e->rad_coef || e->isr_coef) {
            double p1=1+r[5], deltaFactor=p1*p1;
            double fx=delta_px/(k[n-2]*ds),fy=delta_py/(k[n-2]*ds);
            double f2=fx*fx+fy*fy;
            double distance=sqrt(1+r[1]*r[1]+r[3]*r[3])*ds;
            px/=p1;py/=p1;
            r[5]-=e->rad_coef*deltaFactor*f2*distance;
            if(e->isr_coef && e->gaussian)
                r[5]-=e->isr_coef*deltaFactor*pow(f2,.75)*sqrt(distance)*e->gaussian(e->rng);
            if(!(1+r[5]>0)) return 0;
            px*=1+r[5];py*=1+r[5];
        }
        if(r[0]<*xmin) *xmin=r[0];
        if(r[0]>*xmax) *xmax=r[0];
    }
    if(e->rad_coef) {
        double p=e->po*(1+r[5]),beta1=p/sqrt(p*p+1);
        r[4]=beta1*(initial_s/beta0+2*path/(beta0+beta1));
    } else r[4]=initial_s+path;
    *xfinal=r[0];
    return ccb_finite(r);
}

/* Output residuals use the same objectives as ccbend_trajectory_error. */
static int ccb_track(double *r, const CCBend *e, int ideal, double *residual) {
    double b[9],a[9],ints[2][8],gk[2],xmin,xmax,xfinal;
    double *part=r, sign=e->angle<0?-1:1;
    double angle=fabs(e->angle), rho=e->length/angle;
    double yaw=e->yaw*sign*(e->edge_flip?-1:1);
    if(sign<0) {r[0]=-r[0];r[1]=-r[1];r[2]=-r[2];r[3]=-r[3];}
    ccb_coefficients(e,b,a,ideal);
    for(int j=0;j<2;j++) {
        int src=e->edge_flip?1-j:j;
        memcpy(ints[j],e->integrals[src],sizeof(ints[j]));
        if(e->edge_flip) {ints[j][0]*=-1;ints[j][3]*=-1;ints[j][5]*=-1;}
        gk[j]=e->gap*e->fint[src];
    }
    double xi=r[0],xpi=r[1];
    switchRbendPlane(&part,1,angle/2-yaw,0);
    r[0]-=e->dx_offset;
    verticalRbendFringe(&part,1,angle/2-yaw,rho,b[1],2*b[2],gk[0],
                       ints[0],(short)sign,0,(short)e->fringe,(short)e->edge_order);
    CCBend body=*e;
    if(ideal) body.rad_coef=body.isr_coef=0;
    if(!ccb_finite(r) || !ccb_body(r,&body,b,a,&xmin,&xmax,&xfinal)) return 0;
    verticalRbendFringe(&part,1,angle/2+yaw,rho,b[1],2*b[2],gk[1],
                       ints[1],(short)sign,1,(short)e->fringe,(short)e->edge_order);
    if(!ideal) r[0]-=e->x_adjust;
    switchRbendPlane(&part,1,angle/2+yaw,0);
    if(residual) {
        residual[0]=e->fringe?xi-r[0]:xmin+xmax;
        residual[1]=(xpi+r[1])/angle;
        residual[2]=e->fringe?r[0]:xfinal;
    }
    if(sign<0) {r[0]=-r[0];r[1]=-r[1];r[2]=-r[2];r[3]=-r[3];}
    return ccb_finite(r);
}

static double ccb_objective(CCBend *e, const double *p, double *res) {
    double r[6]={0};
    if(fabs(p[0])>=1 || fabs(p[1])>=1) return DBL_MAX;
    e->fse_offset=p[0];e->dx_offset=p[1];
    if(!ccb_track(r,e,1,res)) return DBL_MAX;
    return (e->optimize_dx?fabs(res[0]):0)+(e->optimize_fse?fabs(res[1]):0);
}

/* Bounded Nelder-Mead, same residuals as elegant, local state per element. */
static int ccb_prepare(CCBend *e) {
    double p[3][2]={{0,0},{0,0},{0,0}}, val[3],res[3];
    int active[2],n=0;
    if(e->optimize_fse) active[n++]=0;
    if(e->optimize_dx) active[n++]=1;
    for(int i=0;i<=n;i++) {
        p[i][0]=e->fse_offset;p[i][1]=e->dx_offset;
        if(i) p[i][active[i-1]]+=active[i-1]==0?1e-3:1e-4;
        val[i]=ccb_objective(e,p[i],res);
    }
    int converged=!n;
    for(int iter=0;n && iter<6000;iter++) {
        for(int i=0;i<=n;i++) for(int j=i+1;j<=n;j++) if(val[j]<val[i]) {
            double t=val[i];val[i]=val[j];val[j]=t;
            for(int k=0;k<2;k++){t=p[i][k];p[i][k]=p[j][k];p[j][k]=t;}
        }
        if(val[0]<2e-14*fmax(1,e->length)) {converged=1;break;}
        double c[2]={0},q[2],v;
        for(int k=0;k<2;k++) {
            for(int i=0;i<n;i++) c[k]+=p[i][k]/n;
            q[k]=2*c[k]-p[n][k];
        }
        v=ccb_objective(e,q,res);
        if(v<val[0]) {
            double z[2]={3*c[0]-2*p[n][0],3*c[1]-2*p[n][1]};
            double vz=ccb_objective(e,z,res);
            if(vz<v) {memcpy(q,z,sizeof(q));v=vz;}
        } else if(v>=val[n-1]) {
            double z[2];
            for(int k=0;k<2;k++) z[k]=c[k]+.5*((v<val[n]?q[k]:p[n][k])-c[k]);
            double vz=ccb_objective(e,z,res);
            if(vz<fmin(v,val[n])) {memcpy(q,z,sizeof(q));v=vz;}
            else {
                for(int i=1;i<=n;i++) {
                    for(int k=0;k<2;k++) p[i][k]=.5*(p[i][k]+p[0][k]);
                    val[i]=ccb_objective(e,p[i],res);
                }
                continue;
            }
        }
        memcpy(p[n],q,sizeof(q));val[n]=v;
    }
    if(!converged) return 0;
    ccb_objective(e,p[0],res);
    e->x_adjust=n?res[2]:0;
    double r[6]={0};
    /* Reference corrections are computed with error fields zeroed. */
    CCBend nominal=*e;
    nominal.fse=nominal.fse_dipole=nominal.fse_quad=0;
    nominal.rad_coef=nominal.isr_coef=0;
    memset(r,0,sizeof(r));
    if(!ccb_track(r,&nominal,0,NULL)) return 0;
    for(int i=0;i<4;i++) e->reference[i]=r[i];
    e->reference[4]=r[4]-e->length;
    return 1;
}

static int ccb_at_track(double *r, const CCBend *e) {
    double z2=(1+r[4])*(1+r[4])-r[1]*r[1]-r[3]*r[3];
    if(!(z2>0) || !(1+r[4]>0)) return 0;
    double q[6]={r[0],r[1]/sqrt(z2),r[2],r[3]/sqrt(z2),r[5],r[4]};
    if(!ccb_track(q,e,0,NULL)) return 0;
    if(e->optimize_fse || e->optimize_dx) {
        if(e->correction&2) for(int i=0;i<4;i++) q[i]-=e->reference[i];
        if(e->correction&1) q[4]-=e->reference[4];
    }
    double norm=(1+q[5])/sqrt(1+q[1]*q[1]+q[3]*q[3]);
    r[0]=q[0];r[1]=norm*q[1];r[2]=q[2];r[3]=norm*q[3];
    r[4]=q[5];r[5]=q[4]-e->length;
    return ccb_finite(r);
}
#endif
