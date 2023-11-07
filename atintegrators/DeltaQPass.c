#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"

struct elem 
{
  double Alphax;
  double Alphay;
  double Betax;
  double Betay;
  double Qpx;
  double Qpy;
  double A1;
  double A2;
  double A3;  
  /* optional fields */
  double *R1;
  double *R2;
  double *T1;
  double *T2;
};

void DeltaQPass(double *r_in, int num_particles, double alphax, double alphay,
        double betax, double betay, double qpx, double qpy,
        double a1,double a2, double a3,
        const double *T1, const double *T2,
        const double *R1, const double *R2)
{
     /*
     r_in - 6-by-N matrix of initial conditions reshaped into
     1-d array of 6*N elements
     */
    int i;
    double *rtmp;
    double x,xp,y,yp,dpp;
    double jx,jy;
    double gx,gy;
    double dqx,dqy;
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
            x = rtmp[0];
            xp = rtmp[1]/(1.0+dpp);
            y = rtmp[2];
            yp = rtmp[3]/(1.0+dpp);
            
            jx = 0.5*(gx*x*x+2.0*alphax*x*xp+betax*xp*xp);
            jy = 0.5*(gy*y*y+2.0*alphay*y*yp+betay*yp*yp);
            
            dqx = qpx*dpp+a1*jx+a2*jy;
            dqy = qpy*dpp+a2*jx+a3*jy;
            
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
            
            r_in[i*6] = r11x*x+r12x*xp;
            r_in[i*6+1] = (r21x*x+r22x*xp)*(1+dpp);
            r_in[i*6+2] = r11y*y+r12y*yp;
            r_in[i*6+3] = (r21y*y+r22y*yp)*(1+dpp);
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
    if (!Elem) {
        double alphax, alphay, betax, betay, qpx, qpy, a1, a2, a3;
        double  *R1, *R2, *T1, *T2;
        alphax=atGetDouble(ElemData,"Alphax"); check_error();
        alphay=atGetDouble(ElemData,"Alphay"); check_error();
        betax=atGetDouble(ElemData,"Betax"); check_error();
        betay=atGetDouble(ElemData,"Betay"); check_error();
        qpx=atGetDouble(ElemData,"Qpx"); check_error();
        qpy=atGetDouble(ElemData,"Qpy"); check_error();
        a1=atGetDouble(ElemData,"A1"); check_error();
        a2=atGetDouble(ElemData,"A2"); check_error();
        a3=atGetDouble(ElemData,"A3"); check_error();
        /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Alphax=alphax;
        Elem->Alphay=alphay;
        Elem->Betax=betax;
        Elem->Betay=betay;
        Elem->Qpx=qpx;
        Elem->Qpy=qpy;
        Elem->A1=a1;
        Elem->A2=a2;
        Elem->A3=a3;
        /*optional fields*/
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
    DeltaQPass(r_in, num_particles, Elem->Alphax, Elem->Alphay, 
            Elem->Betax, Elem->Betay, Elem->Qpx, Elem->Qpy,
            Elem->A1, Elem->A2, Elem->A3, Elem->T1, Elem->T2, Elem->R1, Elem->R2);
    return Elem;
}

MODULE_DEF(DeltaQPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double alphax, alphay, betax, betay, qpx, qpy, a1, a2, a3;
        double  *R1, *R2, *T1, *T2;
        alphax=atGetDouble(ElemData,"Alphax"); check_error();
        alphay=atGetDouble(ElemData,"Alphay"); check_error();
        betax=atGetDouble(ElemData,"Betax"); check_error();
        betay=atGetDouble(ElemData,"Betay"); check_error();
        qpx=atGetDouble(ElemData,"Qpx"); check_error();
        qpy=atGetDouble(ElemData,"Qpy"); check_error();
        a1=atGetDouble(ElemData,"A1"); check_error();
        a2=atGetDouble(ElemData,"A2"); check_error();
        a3=atGetDouble(ElemData,"A3"); check_error();
        /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        DeltaQPass(r_in, num_particles, alphax, alphay, 
            betax, betay, qpx, qpy,
            a1, a2, a3, T1, T2, R1, R2);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(9,1);
        mxSetCell(plhs[0],0,mxCreateString("Alphax"));
        mxSetCell(plhs[0],1,mxCreateString("Alphay"));
        mxSetCell(plhs[0],2,mxCreateString("Betax"));
        mxSetCell(plhs[0],3,mxCreateString("Betay"));
        mxSetCell(plhs[0],4,mxCreateString("Qpx"));
        mxSetCell(plhs[0],5,mxCreateString("Qpy"));
        mxSetCell(plhs[0],6,mxCreateString("A1"));
        mxSetCell(plhs[0],7,mxCreateString("A2"));
        mxSetCell(plhs[0],8,mxCreateString("A3"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(4,1);
            mxSetCell(plhs[1],0,mxCreateString("T1"));
            mxSetCell(plhs[1],1,mxCreateString("T2"));
            mxSetCell(plhs[1],2,mxCreateString("R1"));
            mxSetCell(plhs[1],3,mxCreateString("R2"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/

