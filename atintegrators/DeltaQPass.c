#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"

struct elem 
{
  double *Alpha;
  double *Beta;
  double *ChromX;
  double *ChromY;
  int Chrom_MaxOrder;
  double *Detuning; 
  /* optional fields */
  double *Dispersion;
  double *R1;
  double *R2;
  double *T1;
  double *T2;
  double *Alphac;
  int Alphac_MaxOrder;
};

void DeltaQPass(double *r_in, int num_particles, double *alpha,
        double *beta, double *dispersion, double *chromx, double *chromy,
        int chrom_maxorder, double *detuning,
        double *alphac, int alphac_maxorder, double circumference,
        const double *T1, const double *T2,
        const double *R1, const double *R2)
{
     /*
     r_in - 6-by-N matrix of initial conditions reshaped into
     1-d array of 6*N elements
     */
    double alphax = alpha[0];
    double alphay = alpha[1];
    double betax = beta[0];
    double betay = beta[1];
    double dx = 0.0;
    double dy = 0.0;
    double dpx = 0.0;
    double dpy = 0.0;
    if(dispersion){
        dx = dispersion[0];
        dy = dispersion[1];
        dpx = dispersion[2];
        dpy = dispersion[3];
    }
    double a1 = detuning[0];
    double a2 = detuning[1];
    double a3 = detuning[2];

    int i,iq;
    double *rtmp;
    double x,xp,y,yp,dpp, tmpdp;
    double jx,jy;
    double gx,gy;
    double dqx_chrom, dqy_chrom, factorial;
    double dqx,dqy;
    double dct;
    double r11x,r12x,r21x,r22x;
    double r11y,r12y,r21y,r22y;
    bool useT1 = (T1 != NULL);
    bool useT2 = (T2 != NULL);
    bool useR1 = (R1 != NULL);
    bool useR2 = (R2 != NULL);
    double cxy, sxy;
    
    gx = (1+alphax*alphax)/betax;
    gy = (1+alphay*alphay)/betay;
    
    for(i=0; i<num_particles; i++) {
        rtmp = r_in+i*6;
        if(!atIsNaN(rtmp[0])) {
            /*  misalignment at entrance  */
            if (useT1) ATaddvv(rtmp, T1);
            if (useR1) ATmultmv(rtmp, R1);
            dpp = rtmp[4];
            x = rtmp[0]-dx*dpp;
            xp = (rtmp[1]-dpx*dpp)/(1.0+dpp);
            y = rtmp[2]-dy*dpp;
            yp = (rtmp[3]-dpy*dpp)/(1.0+dpp);
            
            jx = 0.5*(gx*x*x+2.0*alphax*x*xp+betax*xp*xp);
            jy = 0.5*(gy*y*y+2.0*alphay*y*yp+betay*yp*yp);
            
            dqx_chrom = 0.0 ; dqy_chrom = 0.0; factorial=1.0; tmpdp = dpp;
            for(iq=0;iq<chrom_maxorder+1; iq++) {
                factorial *= iq + 1;
                dqx_chrom += chromx[iq] * tmpdp / factorial;
                dqy_chrom += chromy[iq] * tmpdp / factorial;
                tmpdp *= dpp;
            }
            
            dct = 0.0;
            if(alphac && alphac_maxorder>0){
                /*Start at second order*/
                tmpdp = dpp*dpp;
                for(iq=1;iq<alphac_maxorder+1; iq++) {
                    dct += alphac[iq] * tmpdp;
                    tmpdp *= dpp;
                }
                dct *= circumference;
            }

            dqx = dqx_chrom + a1*jx + a2*jy;
            dqy = dqy_chrom + a2*jx + a3*jy;
            
            cxy = cos(TWOPI*dqx);
            sxy = sin(TWOPI*dqx);
            r11x = cxy+alphax*sxy;
            r12x = betax*sxy;
            r21x = -gx*sxy;
            r22x = cxy-alphax*sxy;
            
            cxy = cos(TWOPI*dqy);
            sxy = sin(TWOPI*dqy);
            r11y = cxy+alphay*sxy;
            r12y = betay*sxy;
            r21y = -gy*sxy;
            r22y = cxy-alphay*sxy;
            
            r_in[i*6] = r11x*x+r12x*xp+dx*dpp;
            r_in[i*6+1] = (r21x*x+r22x*xp)*(1+dpp)+dpx*dpp;
            r_in[i*6+2] = r11y*y+r12y*yp+dy*dpp;
            r_in[i*6+3] = (r21y*y+r22y*yp)*(1+dpp)+dpy*dpp;
            r_in[i*6+5] += dct;
            /* Misalignment at exit */
            if (useR2) ATmultmv(rtmp, R2);
            if (useT2) ATaddvv(rtmp, T2);
        }    
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                double *r_in, int num_particles, struct parameters *Param)
{
    double circumference = Param->RingLength;
    if (!Elem) {
        int chrom_maxorder, alphac_maxorder;
        double *alpha, *beta, *detuning, *dispersion;
        double  *R1, *R2, *T1, *T2, *chromx, *chromy, *alphac;
        alpha=atGetDoubleArray(ElemData,"Alpha"); check_error();
        beta=atGetDoubleArray(ElemData,"Beta"); check_error();
        chromx=atGetDoubleArray(ElemData,"ChromX"); check_error();
        chromy=atGetDoubleArray(ElemData,"ChromY"); check_error();
        chrom_maxorder=atGetLong(ElemData,"Chrom_MaxOrder"); check_error();
        detuning=atGetDoubleArray(ElemData,"Detuning"); check_error();
        /*optional fields*/
        dispersion=atGetOptionalDoubleArray(ElemData,"Dispersion"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        alphac=atGetOptionalDoubleArray(ElemData,"Alphac"); check_error();
        alphac_maxorder=atGetOptionalLong(ElemData,"Alphac_MaxOrder", 1); check_error();
   
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Alpha=alpha;
        Elem->Beta=beta;
        Elem->Dispersion=dispersion;
        Elem->ChromX=chromx;
        Elem->ChromY=chromy;
        Elem->Chrom_MaxOrder=chrom_maxorder;
        Elem->Detuning=detuning;
        /*optional fields*/
        Elem->Alphac = alphac;
        Elem->Alphac_MaxOrder=alphac_maxorder;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
    DeltaQPass(r_in, num_particles, Elem->Alpha, Elem->Beta,
               Elem->Dispersion, Elem->ChromX, Elem->ChromY,
               Elem->Chrom_MaxOrder, Elem->Detuning,
               Elem->Alphac, Elem->Alphac_MaxOrder, circumference,
               Elem->T1, Elem->T2, Elem->R1, Elem->R2);
    return Elem;
}

MODULE_DEF(DeltaQPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        int chrom_maxorder, alphac_maxorder;
        double *alpha, *beta, *detuning, *dispersion;
        double  *R1, *R2, *T1, *T2, *chromx, *chromy, *alphac;
        alpha=atGetDoubleArray(ElemData,"Alpha"); check_error();
        beta=atGetDoubleArray(ElemData,"Beta"); check_error();
        chromx=atGetDoubleArray(ElemData,"ChromX"); check_error();
        chromy=atGetDoubleArray(ElemData,"ChromY"); check_error();
        chrom_maxorder=atGetLong(ElemData,"Chrom_MaxOrder"); check_error();
        detuning=atGetDoubleArray(ElemData,"Detuning"); check_error();
        /*optional fields*/
        dispersion=atGetOptionalDoubleArray(ElemData,"Dispersion"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        alphac=atGetOptionalDoubleArray(ElemData,"Alphac"); check_error();
        alphac_maxorder=atGetOptionalLong(ElemData,"Alphac_MaxOrder", 1); check_error();
        

      /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        DeltaQPass(r_in, num_particles, alpha, beta, dispersion,
                   chromx, chromy, chrom_maxorder, detuning, alphac,
                   alphac_maxorder, 1.0, T1, T2, R1, R2);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("Alpha"));
        mxSetCell(plhs[0],1,mxCreateString("Beta"));
        mxSetCell(plhs[0],2,mxCreateString("ChromX"));
        mxSetCell(plhs[0],3,mxCreateString("ChromY"));
        mxSetCell(plhs[0],4,mxCreateString("Chrom_MaxOrder"));
        mxSetCell(plhs[0],5,mxCreateString("Detuning"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(7,1);
            mxSetCell(plhs[1],0,mxCreateString("T1"));
            mxSetCell(plhs[1],1,mxCreateString("T2"));
            mxSetCell(plhs[1],2,mxCreateString("R1"));
            mxSetCell(plhs[1],3,mxCreateString("R2"));
            mxSetCell(plhs[1],4,mxCreateString("Alphac"));
            mxSetCell(plhs[1],5,mxCreateString("Alphac_MaxOrder"));
            mxSetCell(plhs[1],6,mxCreateString("Dispersion"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/

