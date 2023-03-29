/* Version modified by Laurent S. Nadolski - Synchrotron SOLEIL
 * April 6th, 2007
 * Modifications
 * 1 - english comments added
 * 2 - debug information display if fourth argument is 1
 * 3 - optional outputs: amplitude and phase
 *
 * known NAFF bugs: data length has to be at least 64
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#ifndef OCTAVE
#include "matrix.h"
#endif
/* #include "gpfunc.h" */
#include "modnaff.h"
#include "complexe.h"
/* #include <sys/ddi.h> */

/* Get ready for R2018a C matrix API */
#ifndef mxGetDoubles
#define mxGetDoubles mxGetPr
#define mxSetDoubles mxSetPr
typedef double mxDouble;
#endif

/* Input Arguments */
#define	Y_IN    prhs[0]
#define	YP_IN   prhs[1]
#define	WIN_IN  prhs[2]
#define	NFREQ_IN  prhs[3]
#define	DEBUG_IN  prhs[4]

/* Output Arguments */
#define	NU_OUT plhs[0]
#define	AMPLITUDE_OUT plhs[1]
#define	PHASE_OUT plhs[2]

unsigned int call_naff(double *ydata, double *ypdata, int ndata, double  *nu_out,
double *amplitude_out, double *phase_out, int win, int nfreq, int debug)
{
    int i, iCpt;
    unsigned int numfreq;
    double *d_in, *d_out;
    t_complexe *c_in;
    
    /* ndata is truncated to be a multiple of 6 if is not yet)
     * ndata-1 to make sure valid data is read when 
     * writing ydata/ypdata to g_NAFVariable.ZTABS */
    ndata = 6*(int )((ndata-1)/6.0);
    
    g_NAFVariable.DTOUR=2*M_PI; /* size of a "cadran" */
    g_NAFVariable.XH=1;         /* step */
    g_NAFVariable.T0=0;         /* time t0 */
    g_NAFVariable.NTERM=nfreq;     /* max term to find */
    g_NAFVariable.KTABS=ndata;  /* number of data : must be a multiple of 6 */
    g_NAFVariable.m_pListFen=NULL; /*no window*/
    g_NAFVariable.TFS=NULL;    /* will contain frequency */
    g_NAFVariable.ZAMP=NULL;   /* will contain amplitude */
    g_NAFVariable.ZTABS=NULL;  /* will contain data to analyze */
    
 /*internal use in naf */
    g_NAFVariable.NERROR=0;
    g_NAFVariable.ICPLX=1;
    g_NAFVariable.IPRT=-1; /*0*/
    g_NAFVariable.NFPRT=stdout; /*NULL;*/
    g_NAFVariable.NFS=0;
    g_NAFVariable.IW=win;
    g_NAFVariable.ISEC=1;
    g_NAFVariable.EPSM=2.2204e-16;
    g_NAFVariable.UNIANG=0;
    g_NAFVariable.FREFON=0;
    g_NAFVariable.ZALP=NULL;
    g_NAFVariable.m_iNbLineToIgnore=1; /*unused*/
    g_NAFVariable.m_dneps=1.E100;
    g_NAFVariable.m_bFSTAB=FALSE; /*unused*/
 /*end of internal use in naf */
    
    naf_initnaf();
    
 /*Transform initial data to complex data since algorithm is optimized for cx data*/
    for (i=0;i<=ndata;i++) {
        g_NAFVariable.ZTABS[i].reel = ydata[i];
        g_NAFVariable.ZTABS[i].imag = ypdata[i];
    }
    
    /*Frequency map analysis*/
    /* look for at most nfreq first fundamental frequencies */
    naf_mftnaf(nfreq,fabs(g_NAFVariable.FREFON)/g_NAFVariable.m_dneps);

    /* number of identified fondamental frequencies */
    numfreq = (unsigned int) g_NAFVariable.NFS;
    
	/* print out results */
    if (debug == 1) {
        mexPrintf("*** NAFF results ***\n");
        mexPrintf("NFS = %d\n",numfreq);
        
        d_in=g_NAFVariable.TFS+1;
        c_in=g_NAFVariable.ZAMP+1;
        for (iCpt=1;iCpt<=numfreq; iCpt++) {
            mexPrintf("AMPL=% 9.6e+i*% 9.6e abs(AMPL)=% 9.6e arg(AMPL)=% 9.6e FREQ=% 9.6e\n",
            c_in->reel,c_in->imag,i_compl_module(*c_in),i_compl_angle(*c_in),*d_in++);
            c_in++;
        }
    }
    
    d_in=g_NAFVariable.TFS+1;
    c_in=g_NAFVariable.ZAMP+1;
    for (iCpt=1;iCpt<=numfreq; iCpt++) {
        *nu_out++ = *d_in++;
        *amplitude_out++ = i_compl_module(*c_in);
        *phase_out++ = i_compl_angle(*c_in);
        c_in++;
    }
    
    
    /*free memory*/
    naf_cleannaf();
    
    /* return number of fundamental frequencies */
    return numfreq;
}

#define min(a, b)       ((a) < (b) ? (a) : (b))
#define max(a, b)       ((a) < (b) ? (b) : (a))
#define NFREQMAX 10 /* maximum number of frequencies to look for */

/*  MATLAB TO C-CALL LINKING FUNCTION  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double   *nu, *amplitude, *phase;
    unsigned int  i, m, n, m2, n2, numfreq;
    int win = 0;
    int debug = 0;          /* No debug */
    int nfreq = NFREQMAX;	/* maximum number of frequencies value per default */
    
    if ((nrhs < 2) || (nrhs >5)) {
        mexErrMsgTxt("Requires 2 to 5 input arguments");
    }
    
    if (nrhs >= 3){ /* windowing */
        m = mxGetM(WIN_IN);
        n = mxGetN(WIN_IN);
        
        if (!mxIsNumeric(WIN_IN) || mxIsComplex(WIN_IN) ||
        mxIsSparse(WIN_IN)  ||  (max(m,n) != 1) || (min(m,n) != 1))
        {
            mexErrMsgTxt("CALCNAFF requires that Window be a scalar.");
        }
        
        win = mxGetScalar(WIN_IN);
    }
    
    if (nrhs >= 4){ /* user fequency number */
        m = mxGetM(NFREQ_IN);
        n = mxGetN(NFREQ_IN);
        
        if (!mxIsNumeric(NFREQ_IN) || mxIsComplex(NFREQ_IN) ||
        mxIsSparse(NFREQ_IN)  ||  (max(m,n) != 1) || (min(m,n) != 1))
        {
            mexErrMsgTxt("CALCNAFF requires that frequency number be a scalar.");
        }

        nfreq = (int )mxGetScalar(NFREQ_IN);
    }
    
    if (nrhs >= 5){ /* debugging flag */
        m = mxGetM(DEBUG_IN);
        n = mxGetN(DEBUG_IN);
        
        if (!mxIsNumeric(DEBUG_IN) || mxIsComplex(DEBUG_IN) ||
        mxIsSparse(DEBUG_IN)  ||  (max(m,n) != 1) || (min(m,n) != 1))
        {
            mexErrMsgTxt("CALCNAFF requires that dubbing flag be a scalar.");
        }
        
        debug = (int )mxGetScalar(DEBUG_IN);
        
    }
        
  /* Check the dimensions of Y.  Y can be >6 X 1 or 1 X >6. */
    
    m = mxGetM(Y_IN);
    n = mxGetN(Y_IN);
    
    if (!mxIsNumeric(Y_IN) || mxIsComplex(Y_IN) ||
    mxIsSparse(Y_IN)  || !mxIsDouble(Y_IN) ||
    (max(m,n) < 66) || (min(m,n) != 1)) {
        mexErrMsgTxt("CALCNAFF requires that Y be a >= 66 x 1 vector.");
    }
    
    /* Check the dimensions of YP.  YP must have same size as Y */    
    m2 = mxGetM(YP_IN);
    n2 = mxGetN(YP_IN);
    
    if (!mxIsNumeric(YP_IN) || mxIsComplex(YP_IN) ||
    mxIsSparse(YP_IN)  || !mxIsDouble(YP_IN) ||
    (m2 != m) || (n2 != n)) {
        mexErrMsgTxt("CALCNAFF requires that YP has the same size as Y.");
    }
    
    /* Dynamic memory allocation */
    nu = (double *) mxMalloc(nfreq*sizeof(double));
    amplitude = (double *) mxMalloc(nfreq*sizeof(double));
    phase = (double *) mxMalloc(nfreq*sizeof(double));
    
    /* call subroutine that calls routine for all NAFF computation */
    numfreq = call_naff(mxGetDoubles(Y_IN),mxGetDoubles(YP_IN),(int )max(m,n),
            nu, amplitude, phase,  win, nfreq, debug);
    
    NU_OUT = mxCreateDoubleMatrix(numfreq, 1, mxREAL);
    memcpy(mxGetDoubles(NU_OUT), nu, numfreq*sizeof(double));
    
    if (nlhs >= 2){ /* amplitudes */
        AMPLITUDE_OUT = mxCreateDoubleMatrix(numfreq, 1, mxREAL);
        memcpy(mxGetDoubles(AMPLITUDE_OUT), amplitude, numfreq*sizeof(double));
    }
    
    if (nlhs >= 3){ /* phases */
        PHASE_OUT = mxCreateDoubleMatrix(numfreq, 1, mxREAL);
        memcpy(mxGetDoubles(PHASE_OUT), phase, numfreq*sizeof(double));
    }
    
    /* free dynamically memory allocation */
    mxFree(nu);
    mxFree(amplitude);
    mxFree(phase);
    
    return;
} /* end of mexFunction */
