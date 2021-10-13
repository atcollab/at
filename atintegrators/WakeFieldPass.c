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
  double fx;
  double fy;       
  double fz;
  double fqx;
  double fqy;
  double circumference;
  double *waketableT;
  double *waketableDX;
  double *waketableDY;
  double *waketableQX;
  double *waketableQY;
  double *waketableZ;
  double *turnhistoryX;
  double *turnhistoryY;
  double *turnhistoryZ;
  double *turnhistoryW;
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
    double fx = Elem->fx;
    double fy = Elem->fy;
    double fqx = Elem->fqx;
    double fqy = Elem->fqy;
    double fz = Elem->fz;
    double circumference = Elem->circumference;
    double *waketableT = Elem->waketableT;
    double *waketableDX = Elem->waketableDX;
    double *waketableDY = Elem->waketableDY;
    double *waketableQX = Elem->waketableQX;
    double *waketableQY = Elem->waketableQY;
    double *waketableZ = Elem->waketableZ;
    double *turnhistoryX = Elem->turnhistoryX;
    double *turnhistoryY = Elem->turnhistoryY;
    double *turnhistoryZ = Elem->turnhistoryZ;
    double *turnhistoryW = Elem->turnhistoryW;
    double *z_cuts = Elem->z_cuts;    

    size_t sz = 5*nslice*sizeof(double) + num_particles*sizeof(int);
    int i,ii;
    double *rtmp;

    int *pslice;
    double *kx;
    double *ky;
    double *kx2;
    double *ky2;
    double *kz;

    double wi, dx, dy, ds;
    int index;

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

    rotate_table_history(nturns,nslice,turnhistoryX,turnhistoryY,turnhistoryZ,turnhistoryW,circumference);
    slice_bunch(r_in,num_particles,nslice,nturns,turnhistoryX,turnhistoryY,turnhistoryZ,turnhistoryW,pslice,z_cuts);

    for(i=nslice*(nturns-1);i<nslice*nturns;i++){        
        if(turnhistoryW[i]>0.0){
            for (ii=0;ii<nslice*nturns;ii++){
                ds = turnhistoryZ[ii]-turnhistoryZ[i];
                if(turnhistoryW[ii]>0.0 && -ds>=waketableT[0] && -ds<waketableT[nelem-1]){
                    wi = turnhistoryW[ii];
                    dx = turnhistoryX[ii];
                    dy = turnhistoryY[ii];
                    index = binarySearch(waketableT,-ds,nelem,0,0);              
                    if(waketableDX)kx[i-nslice*(nturns-1)] += dx*fx*wi*getTableWake(waketableDX,waketableT,-ds,index);
                    if(waketableDY)ky[i-nslice*(nturns-1)] += dy*fy*wi*getTableWake(waketableDY,waketableT,-ds,index);
                    if(waketableQX)kx2[i-nslice*(nturns-1)] += fqx*wi*getTableWake(waketableQX,waketableT,-ds,index);
                    if(waketableQY)ky2[i-nslice*(nturns-1)] += fqy*wi*getTableWake(waketableQY,waketableT,-ds,index);
                    if(waketableZ) kz[i-nslice*(nturns-1)] += fz*wi*getTableWake(waketableZ,waketableT,-ds,index);
                }            
            }
        }
    }
    
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
        double normfactx, normfacty, normfactz;
        double *waketableT;
        double *waketableDX;
        double *waketableDY;
        double *waketableQX;
        double *waketableQY;
        double *waketableZ;
        double *turnhistoryX;
        double *turnhistoryY;
        double *turnhistoryZ;
        double *turnhistoryW;
        double *z_cuts;

        nslice=atGetLong(ElemData,"Nslice"); check_error();
        nelem=atGetLong(ElemData,"Nelem"); check_error();  
        nturns=atGetLong(ElemData,"Nturns"); check_error();      
        intensity=atGetDouble(ElemData,"Intensity"); check_error();
        circumference=atGetDouble(ElemData,"Circumference"); check_error();
        wakefact=atGetDouble(ElemData,"Wakefact"); check_error();
        waketableT=atGetDoubleArray(ElemData,"WakeT"); check_error();
        turnhistoryX=atGetDoubleArray(ElemData,"TurnHistoryX"); check_error();
        turnhistoryY=atGetDoubleArray(ElemData,"TurnHistoryY"); check_error();
        turnhistoryZ=atGetDoubleArray(ElemData,"TurnHistoryZ"); check_error();
        turnhistoryW=atGetDoubleArray(ElemData,"TurnHistoryW"); check_error();
        /*optional attributes*/
        normfactx=atGetOptionalDouble(ElemData,"Normfactx",1.0); check_error();
        normfacty=atGetOptionalDouble(ElemData,"Normfacty",1.0); check_error();
        normfactz=atGetOptionalDouble(ElemData,"Normfactz",1.0); check_error();
        waketableDX=atGetOptionalDoubleArray(ElemData,"WakeDX"); check_error();
        waketableDY=atGetOptionalDoubleArray(ElemData,"WakeDY"); check_error();
        waketableQX=atGetOptionalDoubleArray(ElemData,"WakeQX"); check_error();
        waketableQY=atGetOptionalDoubleArray(ElemData,"WakeQY"); check_error();
        waketableZ=atGetOptionalDoubleArray(ElemData,"WakeZ"); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"z_cuts"); check_error();

        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->nslice=nslice;
        Elem->nelem=nelem;
        Elem->nturns=nturns;
        Elem->circumference=circumference;
        Elem->fx=intensity*wakefact*normfactx;
        Elem->fy=intensity*wakefact*normfacty;
        Elem->fqx=intensity*wakefact*normfactx;
        Elem->fqy=intensity*wakefact*normfacty;
        Elem->fz=intensity*wakefact*normfactz;
        Elem->waketableT=waketableT;
        Elem->waketableDX=waketableDX;
        Elem->waketableDY=waketableDY;
        Elem->waketableQX=waketableQX;
        Elem->waketableQY=waketableQY;
        Elem->waketableZ=waketableZ;
        Elem->turnhistoryX=turnhistoryX;
        Elem->turnhistoryY=turnhistoryY;
        Elem->turnhistoryZ=turnhistoryZ;
        Elem->turnhistoryW=turnhistoryW;
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
        double normfactx, normfacty, normfactz;
        double *waketableT;
        double *waketableDX;
        double *waketableDY;
        double *waketableQX;
        double *waketableQY;
        double *waketableZ;
        double *turnhistoryX;
        double *turnhistoryY;
        double *turnhistoryZ;
        double *turnhistoryW;
        double *z_cuts;
        
        nslice=atGetLong(ElemData,"Nslice"); check_error();
        nelem=atGetLong(ElemData,"Nelem"); check_error();
        nturns=atGetLong(ElemData,"Nturns"); check_error();
        intensity=atGetDouble(ElemData,"Intensity"); check_error();
        circumference=atGetDouble(ElemData,"Circumference"); check_error();
        wakefact=atGetDouble(ElemData,"Wakefact"); check_error();
        waketableT=atGetDoubleArray(ElemData,"WakeT"); check_error();
        turnhistoryX=atGetDoubleArray(ElemData,"TurnHistoryX"); check_error();
        turnhistoryY=atGetDoubleArray(ElemData,"TurnHistoryY"); check_error();
        turnhistoryZ=atGetDoubleArray(ElemData,"TurnHistoryZ"); check_error();
        turnhistoryW=atGetDoubleArray(ElemData,"TurnHistoryW"); check_error();
        /*optional attributes*/
        normfactx=atGetOptionalDouble(ElemData,"Normfactx",1.0); check_error();
        normfacty=atGetOptionalDouble(ElemData,"Normfacty",1.0); check_error();
        normfactz=atGetOptionalDouble(ElemData,"Normfactz",1.0); check_error();
        waketableDX=atGetOptionalDoubleArray(ElemData,"WakeDX"); check_error();
        waketableDY=atGetOptionalDoubleArray(ElemData,"WakeDY"); check_error();
        waketableQX=atGetOptionalDoubleArray(ElemData,"WakeQX"); check_error();
        waketableQY=atGetOptionalDoubleArray(ElemData,"WakeQY"); check_error();
        waketableZ=atGetOptionalDoubleArray(ElemData,"WakeZ"); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"z_cuts"); check_error();
        
        Elem->nslice=nslice;
        Elem->nelem=nelem;
        Elem->nturns=nturns;
        Elem->circumference=circumference;
        Elem->fx=intensity*wakefact*normfactx;
        Elem->fy=intensity*wakefact*normfacty;
        Elem->fqx=intensity*wakefact*normfactx;
        Elem->fqy=intensity*wakefact*normfacty;
        Elem->fz=intensity*wakefact*normafactz;
        Elem->waketableT=waketableT;
        Elem->waketableDX=waketableDX;
        Elem->waketableDY=waketableDY;
        Elem->waketableQX=waketableQX;
        Elem->waketableQY=waketableQY;
        Elem->waketableZ=waketableZ;
        Elem->turnhistoryX=turnhistoryX;
        Elem->turnhistoryY=turnhistoryY;
        Elem->turnhistoryZ=turnhistoryZ;
        Elem->turnhistoryW=turnhistoryW;
        Elem->z_cuts=z_cuts;

        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix: particle array");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        WakeFieldPass(r_in, num_particles, Elem);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(11,1);
        mxSetCell(plhs[0],0,mxCreateString("Nelem"));
        mxSetCell(plhs[0],1,mxCreateString("Nslice"));
        mxSetCell(plhs[0],2,mxCreateString("Nturns"));
        mxSetCell(plhs[0],3,mxCreateString("Intensity"));
        mxSetCell(plhs[0],4,mxCreateString("Circumference"));
        mxSetCell(plhs[0],5,mxCreateString("Wakefact"));
        mxSetCell(plhs[0],6,mxCreateString("WakeT"));
        mxSetCell(plhs[0],7,mxCreateString("TurnHistoryX"));
        mxSetCell(plhs[0],8,mxCreateString("TurnHistoryY"));
        mxSetCell(plhs[0],9,mxCreateString("TurnHistoryZ"));
        mxSetCell(plhs[0],10,mxCreateString("TurnHistoryW"));

        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(9,1); /* No optional fields */
            mxSetCell(plhs[0],0,mxCreateString("WakeDX"));
            mxSetCell(plhs[0],1,mxCreateString("WakeDY"));
            mxSetCell(plhs[0],2,mxCreateString("WakeQX"));
            mxSetCell(plhs[0],3,mxCreateString("WakeQY"));
            mxSetCell(plhs[0],4,mxCreateString("WakeZ"));
            mxSetCell(plhs[0],5,mxCreateString("Normfactx"));
            mxSetCell(plhs[0],6,mxCreateString("Normfacty"));
            mxSetCell(plhs[0],7,mxCreateString("Normfactz"));
            mxSetCell(plhs[0],8,mxCreateString("z_cuts"));
        }
    }else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif
