/* IdTablePass.c
 * Accelerator Toolbox
 * Created: 13/11/08
 * Z.Mart?? zeus@cells.es
 *
 * Based in the matlab routine:
 * WigTablePass.m - The tracking table is described in
 * P. Elleaume, "A new approach to the electron beam dynamics in undulators
 * and wigglers", EPAC92.
 *
 */

#include <math.h>
#include "mex.h"
#include "elempass.h"
#include "atlalib.c"
#include "interpolate.c"
#include "matrix.h"

double *GLOBAL_x,*GLOBAL_y,*GLOBAL_xkick1,*GLOBAL_ykick1,*GLOBAL_xkick,*GLOBAL_ykick,*GLOBAL_xkick2,*GLOBAL_ykick2;
int GLOBAL_m,GLOBAL_n;

/*Definition of the interpolated functions*/
static double Map_x(double x,double y)
{
    double f;
    /*cubic interpolation*/
    /*splin2(GLOBAL_y,GLOBAL_x,GLOBAL_xkick,GLOBAL_xkick2,GLOBAL_n,GLOBAL_m,y,x,&f);*/
    
    /*biliniar interpolation*/
    linint(GLOBAL_y,GLOBAL_x,GLOBAL_xkick,GLOBAL_m,GLOBAL_n,y,x,&f);
    return f;
}

static double Map_y(double x,double y)
{
    double f;
    /*cubic interpolation*/
    /*splin2(GLOBAL_y,GLOBAL_x,GLOBAL_ykick,GLOBAL_ykick2,GLOBAL_m,GLOBAL_n,y,x,&f);*/
    
    /*biliniar interpolation*/
    linint(GLOBAL_y,GLOBAL_x,GLOBAL_ykick,GLOBAL_m,GLOBAL_n,y,x,&f);
    return f;
}

static double Map1_x(double x,double y)
{
    double f;
    /*cubic interpolation*/
    /*splin2(GLOBAL_y,GLOBAL_x,GLOBAL_xkick1,GLOBAL_xkick2,GLOBAL_n,GLOBAL_m,y,x,&f);*/
    
    /*biliniar interpolation*/
    linint(GLOBAL_y,GLOBAL_x,GLOBAL_xkick1,GLOBAL_m,GLOBAL_n,y,x,&f);
    return f;
}

static double Map1_y(double x,double y)
{
    double f;
    /*cubic interpolation*/
    /*splin2(GLOBAL_y,GLOBAL_x,GLOBAL_ykick1,GLOBAL_ykick2,GLOBAL_m,GLOBAL_n,y,x,&f);*/
    
    /*biliniar interpolation*/
    linint(GLOBAL_y,GLOBAL_x,GLOBAL_ykick1,GLOBAL_m,GLOBAL_n,y,x,&f);
    return f;
}
/*
static void markaslost(double *r6,int idx)
{
    r6[idx] = mxGetInf();
}
*/
/* Set T1, T2, R1, R2 to NULL pointers to ignore misalignmets*/
void IdKickMapModelPass(double *r, double le, double *xkick1, double *ykick1, double *xkick, double *ykick, double *x, double *y,int n,int m, int Nslice, double *T1, double *T2, double *R1, double *R2, int num_particles)
{
    double *r6,deltaxp,deltayp,deltaxp1,deltayp1,*limitsptr;
    int c;
    bool usexkick1 = (xkick1 != NULL);
    bool useykick1 = (ykick1 != NULL);
    bool useT1 = (T1 != NULL);
    bool useT2 = (T2 != NULL);
    bool useR1 = (R1 != NULL);
    bool useR2 = (R2 != NULL);
    double L1 = le/(2*Nslice);
    
    /*Act as AperturePass*/
    limitsptr=(double*)mxCalloc(4,sizeof(double));
    limitsptr[0]=x[0];
    limitsptr[1]=x[n-1];
    limitsptr[2]=y[0];
    limitsptr[3]=y[m-1];
    
    /*globalize*/
    
    /* For cubic interpolation only*/
    
    /*GLOBAL_xkick2=(double*)mxCalloc(n*m,sizeof(double));
     * GLOBAL_ykick2=(double*)mxCalloc(n*m,sizeof(double));
     * splie2(y,x,xkick,m,n,GLOBAL_xkick2);
     * splie2(y,x,ykick,m,n,GLOBAL_ykick2); */
    
    GLOBAL_x=x;
    GLOBAL_y=y;
    GLOBAL_xkick1=xkick1;
    GLOBAL_ykick1=ykick1;
    GLOBAL_xkick=xkick;
    GLOBAL_ykick=ykick;
    GLOBAL_m=m; /* y used as colums*/
    GLOBAL_n=n; /* x used as rows*/
    
    for (c=0; c<num_particles; c++) {
        r6 = r+c*6;
        
        if(!mxIsNaN(r6[0]) & mxIsFinite(r6[4])) {
            /*
             * function bend6 internally calculates the square root
             * of the energy deviation of the particle
             * To protect against DOMAIN and OVERFLOW error, check if the
             * fifth component of the phase spacevector r6[4] is finite
             */
            if (r6[0]<limitsptr[0] || r6[0]>limitsptr[1])
                markaslost(r6,0);
            else if (r6[2]<limitsptr[2] || r6[2]>limitsptr[3])
                markaslost(r6,2);
            else {
                /* Misalignment at entrance */
                if (useT1) ATaddvv(r6,T1);
                if (useR1) ATmultmv(r6,R1);
                /*Tracking in the main body*/
                for (m=0; m<Nslice; m++) { /* Loop over slices*/
                    ATdrift6(r6,L1);
                    if (!mxIsNaN(r6[0])&&!mxIsNaN(r6[2])) {
                        /*The kick from IDs varies quadratically, not linearly, with energy.   */
                        deltaxp = (1.0/Nslice)*Map_x(r6[0],r6[2])/(1.0+r6[4]);
                        deltayp = (1.0/Nslice)*Map_y(r6[0],r6[2])/(1.0+r6[4]);
                        if(usexkick1)  deltaxp1 = (1.0/Nslice)*Map1_x(r6[0],r6[2]);
                        if(useykick1)  deltayp1 = (1.0/Nslice)*Map1_y(r6[0],r6[2]);
                        r6[1] = r6[1] + deltaxp;
                        if(usexkick1) r6[1]=r6[1]+deltaxp1;
                        r6[3] = r6[3] + deltayp;
                        if(useykick1) r6[3]= r6[3] + deltayp1;
                    }
                    ATdrift6(r6,L1);
                }
                /* Misalignment at exit */
                if (useR2) ATmultmv(r6,R2);
                if (useT2) ATaddvv(r6,T2);
            }
        }
    }
}

#ifndef NOMEX

#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
        double *r_in,int num_particles,int mode)
#define NUM_FIELDS_2_REMEMBER 12
{
    int n_map,m_map;
    int Nslice;
    double le;
    double *pr1, *pr2, *pt1, *pt2, *xkick, *ykick, *xkick1, *ykick1, *x, *y;
    
    switch(mode) {
        case MAKE_LOCAL_COPY:
            /* Find field numbers first
             * Save a list of field number in an array
             * and make returnptr point to that array
             */
            FieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "Length");
            FieldNumbers[1] = GetRequiredFieldNumber(ElemData, "xkick");
            FieldNumbers[2] = GetRequiredFieldNumber(ElemData, "ykick");
            FieldNumbers[3] = GetRequiredFieldNumber(ElemData, "xtable");
            FieldNumbers[4] = GetRequiredFieldNumber(ElemData, "ytable");
            FieldNumbers[5] = GetRequiredFieldNumber(ElemData, "Nslice");
            
            /* Optional fields */
            
            FieldNumbers[6] = mxGetFieldNumber(ElemData, "xkick1");
            FieldNumbers[7] = mxGetFieldNumber(ElemData, "ykick1");
            FieldNumbers[8] = mxGetFieldNumber(ElemData, "R1");
            FieldNumbers[9] = mxGetFieldNumber(ElemData, "R2");
            FieldNumbers[10] = mxGetFieldNumber(ElemData, "T1");
            FieldNumbers[11] = mxGetFieldNumber(ElemData, "T2");
            /* Fall through next section... */
            
        case	USE_LOCAL_COPY:
            /* Get fields from MATLAB using field numbers
             * The second argument ponter to the array of field
             * numbers is previously created with
             * BendLinearPass( ..., MAKE_LOCAL_COPY)
             */
            le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
            xkick = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
            ykick = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
            x = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
            y = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
            Nslice = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
            
            n_map = mxGetN(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
            m_map = mxGetM(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
            
            /* Optional fields */
            
            xkick1 = (FieldNumbers[6] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[6])) : NULL;
            ykick1 = (FieldNumbers[7] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[7])) : NULL;
            pr1 = (FieldNumbers[8] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[8])) : NULL;
            pr2 = (FieldNumbers[9] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[9])) : NULL;
            pt1 = (FieldNumbers[10] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[10])) : NULL;
            pt2 = (FieldNumbers[11] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[11])) : NULL;
            break;
    }
    
    IdKickMapModelPass(r_in, le,xkick1,ykick1,xkick,ykick,x,y,n_map,m_map,Nslice,
            pt1, pt2, pr1, pr2, num_particles);
    
    return FieldNumbers;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        double *xkick1, *ykick1, *pr1, *pr2, *pt1, *pt2;
        mxArray *tmpmxptr = GetRequiredField(prhs[0], "xkick");
        
        double le = mxGetScalar(GetRequiredField(prhs[0], "Length"));
        double *xkick =  mxGetPr(tmpmxptr);
        int n_map = mxGetN(tmpmxptr);
        int m_map = mxGetM(tmpmxptr);
        double *ykick =  mxGetPr(GetRequiredField(prhs[0], "ykick"));
        double *x =  mxGetPr(GetRequiredField(prhs[0], "xtable"));
        double *y =  mxGetPr(GetRequiredField(prhs[0], "ytable"));
        int Nslice = (int)mxGetScalar(GetRequiredField(prhs[0], "Nslice"));
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        
        /*Optional fields*/
        
        tmpmxptr = mxGetField(prhs[0],0,"xkick1");
        xkick1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"ykick1");
        ykick1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"R1");
        pr1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"R2");
        pr2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T1");
        pt1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T2");
        pt2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        IdKickMapModelPass(r_in, le, xkick1,ykick1, xkick,ykick,x,y,n_map,m_map,Nslice,
                pt1, pt2, pr1, pr2, num_particles);
    }
    else if (nrhs == 0) {
        /* return list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("xkick"));
        mxSetCell(plhs[0],2,mxCreateString("ykick"));
        mxSetCell(plhs[0],3,mxCreateString("xtable"));
        mxSetCell(plhs[0],4,mxCreateString("ytable"));
        mxSetCell(plhs[0],5,mxCreateString("Nslice"));
        
        if (nlhs>1) {
            /* Required and optional fields */
            plhs[1] = mxCreateCellMatrix(6,1);
            mxSetCell(plhs[1],0,mxCreateString("xkick1"));
            mxSetCell(plhs[1],1,mxCreateString("ykick1"));
            mxSetCell(plhs[1],2,mxCreateString("T1"));
            mxSetCell(plhs[1],3,mxCreateString("T2"));
            mxSetCell(plhs[1],4,mxCreateString("R1"));
            mxSetCell(plhs[1],5,mxCreateString("R2"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif
