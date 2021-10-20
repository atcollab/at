#include "atelem.c"
#include "atimplib.c"
#include <math.h>
#include <float.h>
/*
 * Wake pass method by Simon White.  
 * User may contact simon.white@esrf.fr for questions and comments.
 */

struct elem
{
  int nslice;
  int nelem;
  int nturns;
  double circumference;
  double *normfact;
  double *waketableT;
  double *waketableDX;
  double *waketableDY;
  double *waketableQX;
  double *waketableQY;
  double *waketableZ;
  double *turnhistory;
  double *z_cuts;
};


void WakeFieldPass(double *r_in,int num_particles, struct elem *Elem){   
    /*
     * r_in - 6-by-N matrix of initial conditions reshaped into
     * 1-d array of 6*N elements
     */   
    long nslice = Elem->nslice;
    long nelem = Elem->nelem;
    long nturns = Elem->nturns;
    double circumference = Elem->circumference;
    double *normfact = Elem->normfact;
    double *waketableT = Elem->waketableT;
    double *waketableDX = Elem->waketableDX;
    double *waketableDY = Elem->waketableDY;
    double *waketableQX = Elem->waketableQX;
    double *waketableQY = Elem->waketableQY;
    double *waketableZ = Elem->waketableZ;
    double *turnhistory = Elem->turnhistory;
    double *z_cuts = Elem->z_cuts;    

    size_t sz = 5*nslice*sizeof(double) + num_particles*sizeof(int);
    int i;
    double *rtmp;

    int *pslice;
    double *kx;
    double *ky;
    double *kx2;
    double *ky2;
    double *kz;

    void *buffer = atMalloc(sz);
    double *dptr = (double *) buffer;
    int *iptr;

    kx = dptr; dptr += nslice;
    ky = dptr; dptr += nslice;
    kx2 = dptr; dptr += nslice;
    ky2 = dptr; dptr += nslice;
    kz = dptr; dptr += nslice;

    iptr = (int *) dptr;
    pslice = iptr; iptr += num_particles;

    for (i=0;i<nslice;i++) {
        kx[i]=0.0;
        ky[i]=0.0;
        kx2[i]=0.0;
        ky2[i]=0.0;
        kz[i]=0.0;
    }

    rotate_table_history(nturns,nslice,turnhistory,circumference);
    slice_bunch(r_in,num_particles,nslice,nturns,turnhistory,pslice,z_cuts);
    compute_kicks(nslice,nturns,nelem,turnhistory,waketableT,waketableDX,
                  waketableDY,waketableQX,waketableQY,waketableZ,
                  normfact,kx,ky,kx2,ky2,kz);
    
    for (i=0; i<num_particles; i++) {
        rtmp = r_in+i*6;
        if (!atIsNaN(rtmp[0])) {
            rtmp[4] += kz[pslice[i]];
            rtmp[1] += (kx[pslice[i]]+rtmp[0]*kx2[pslice[i]])*(1+rtmp[4]);
            rtmp[3] += (ky[pslice[i]]+rtmp[2]*ky2[pslice[i]])*(1+rtmp[4]);
        }
    }    
    atFree(buffer);
};


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        long nslice,nelem, nturns;
        double intensity, wakefact, circumference;
        double *normfact;
        double *waketableT;
        double *waketableDX;
        double *waketableDY;
        double *waketableQX;
        double *waketableQY;
        double *waketableZ;
        double *turnhistory;
        double *z_cuts;

        nslice=atGetLong(ElemData,"Nslice"); check_error();
        nelem=atGetLong(ElemData,"Nelem"); check_error();  
        nturns=atGetLong(ElemData,"Nturns"); check_error();      
        intensity=atGetDouble(ElemData,"Intensity"); check_error();
        circumference=atGetDouble(ElemData,"Circumference"); check_error();
        wakefact=atGetDouble(ElemData,"Wakefact"); check_error();
        waketableT=atGetDoubleArray(ElemData,"WakeT"); check_error();
        turnhistory=atGetDoubleArray(ElemData,"TurnHistory"); check_error();
        normfact=atGetDoubleArray(ElemData,"NormFact"); check_error();
        /*optional attributes*/
        waketableDX=atGetOptionalDoubleArray(ElemData,"WakeDX"); check_error();
        waketableDY=atGetOptionalDoubleArray(ElemData,"WakeDY"); check_error();
        waketableQX=atGetOptionalDoubleArray(ElemData,"WakeQX"); check_error();
        waketableQY=atGetOptionalDoubleArray(ElemData,"WakeQY"); check_error();
        waketableZ=atGetOptionalDoubleArray(ElemData,"WakeZ"); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();

        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->nslice=nslice;
        Elem->nelem=nelem;
        Elem->nturns=nturns;
        Elem->circumference=circumference;
        for(int i=0;i<3;i++){
           normfact[i]*=intensity*wakefact;
        }
        Elem->normfact=normfact;
        Elem->waketableT=waketableT;
        Elem->waketableDX=waketableDX;
        Elem->waketableDY=waketableDY;
        Elem->waketableQX=waketableQX;
        Elem->waketableQY=waketableQY;
        Elem->waketableZ=waketableZ;
        Elem->turnhistory=turnhistory;
        Elem->z_cuts=z_cuts;
    }
    WakeFieldPass(r_in,num_particles,Elem);
    return Elem;
}

MODULE_DEF(WakeFieldPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        struct elem El, *Elem=&El;

        long nslice,nelem,nturns;
        double intensity, wakefact, circumference;
        double *normfact;
        double *waketableT;
        double *waketableDX;
        double *waketableDY;
        double *waketableQX;
        double *waketableQY;
        double *waketableZ;
        double *turnhistory;
        double *z_cuts;
        
        nslice=atGetLong(ElemData,"Nslice"); check_error();
        nelem=atGetLong(ElemData,"Nelem"); check_error();
        nturns=atGetLong(ElemData,"Nturns"); check_error();
        intensity=atGetDouble(ElemData,"Intensity"); check_error();
        circumference=atGetDouble(ElemData,"Circumference"); check_error();
        wakefact=atGetDouble(ElemData,"Wakefact"); check_error();
        waketableT=atGetDoubleArray(ElemData,"WakeT"); check_error();
        turnhistory=atGetDoubleArray(ElemData,"TurnHistory"); check_error();
        normfact=atGetDoubleArray(ElemData,"NormFact"); check_error();
        /*optional attributes*/
        waketableDX=atGetOptionalDoubleArray(ElemData,"WakeDX"); check_error();
        waketableDY=atGetOptionalDoubleArray(ElemData,"WakeDY"); check_error();
        waketableQX=atGetOptionalDoubleArray(ElemData,"WakeQX"); check_error();
        waketableQY=atGetOptionalDoubleArray(ElemData,"WakeQY"); check_error();
        waketableZ=atGetOptionalDoubleArray(ElemData,"WakeZ"); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
        
        Elem->nslice=nslice;
        Elem->nelem=nelem;
        Elem->nturns=nturns;
        Elem->circumference=circumference;
        for(int i=0;i<3;i++){
           normfact[i]*=intensity*wakefact;
        }
        Elem->normfact=normfact;
        Elem->waketableT=waketableT;
        Elem->waketableDX=waketableDX;
        Elem->waketableDY=waketableDY;
        Elem->waketableQX=waketableQX;
        Elem->waketableQY=waketableQY;
        Elem->waketableZ=waketableZ;
        Elem->turnhistoryX=turnhistory;
        Elem->z_cuts=z_cuts;

        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix: particle array");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        WakeFieldPass(r_in, num_particles, Elem);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(9,1);
        mxSetCell(plhs[0],0,mxCreateString("Nelem"));
        mxSetCell(plhs[0],1,mxCreateString("Nslice"));
        mxSetCell(plhs[0],2,mxCreateString("Nturns"));
        mxSetCell(plhs[0],3,mxCreateString("Intensity"));
        mxSetCell(plhs[0],4,mxCreateString("Circumference"));
        mxSetCell(plhs[0],5,mxCreateString("Wakefact"));
        mxSetCell(plhs[0],6,mxCreateString("WakeT"));
        mxSetCell(plhs[0],7,mxCreateString("TurnHistory"));
        mxSetCell(plhs[0],8,mxCreateString("Normfact"));

        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(6,1); /* No optional fields */
            mxSetCell(plhs[0],0,mxCreateString("WakeDX"));
            mxSetCell(plhs[0],1,mxCreateString("WakeDY"));
            mxSetCell(plhs[0],2,mxCreateString("WakeQX"));
            mxSetCell(plhs[0],3,mxCreateString("WakeQY"));
            mxSetCell(plhs[0],4,mxCreateString("WakeZ"));
            mxSetCell(plhs[0],5,mxCreateString("ZCuts"));
        }
    }else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif
