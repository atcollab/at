#include "atelem.c"
#include "atimplib.c"
#include <math.h>
#include <float.h>
/*
 * Impedance pass method by Simon White.  
 * User may contact simon.white@esrf.fr for questions and comments.
 */

struct elem
{
  int nslice;
  int nelem;
  double fx;
  double fy;       
  double fz;
  double fqx;
  double fqy;
  double *waketableT;
  double *waketableDX;
  double *waketableDY;
  double *waketableQX;
  double *waketableQY;
  double *waketableZ;
};


void impedance_tablePass(double *r_in,int num_particles, struct elem *Elem){   
    /*
     * r_in - 6-by-N matrix of initial conditions reshaped into
     * 1-d array of 6*N elements
     */   
    long nslice = Elem->nslice;
    long nelem = Elem->nelem;
    double fx = Elem->fx;
    double fy = Elem->fy;
    double fqx = Elem->fqx;
    double fqy = Elem->fqy;
    double fz = Elem->fz;
    double *waketableT = Elem->waketableT;
    double *waketableDX = Elem->waketableDX;
    double *waketableDY = Elem->waketableDY;
    double *waketableQX = Elem->waketableQX;
    double *waketableQY = Elem->waketableQY;
    double *waketableZ = Elem->waketableZ;
    
    size_t sz = 9*nslice*sizeof(double) + (nslice+num_particles)*sizeof(int);
    int i,ii;
    double *rtmp;

    double *weight;
    double *xpos;
    double *ypos;
    double *zpos;
    int *countslc;
    int *pslice;
    double *kx;
    double *ky;
    double *kx2;
    double *ky2;
    double *kz;

    void *buffer = atMalloc(sz);
    double *dptr = (double *) buffer;
    int *iptr;

    weight = dptr; dptr += nslice;
    xpos = dptr; dptr += nslice;
    ypos = dptr; dptr += nslice;
    zpos = dptr; dptr += nslice;
    kx = dptr; dptr += nslice;
    ky = dptr; dptr += nslice;
    kx2 = dptr; dptr += nslice;
    ky2 = dptr; dptr += nslice;
    kz = dptr; dptr += nslice;

    iptr = (int *) dptr;
    countslc = iptr; iptr += nslice;
    pslice = iptr; iptr += num_particles;

    for (i=0;i<nslice;i++) {
        countslc[i]=0;
        xpos[i]=0.0;
        ypos[i]=0.0;
        zpos[i]=0.0;
        kx[i]=0.0;
        ky[i]=0.0;
        kx2[i]=0.0;
        ky2[i]=0.0;
        kz[i]=0.0;
    }

    double *bounds;
    bounds = getbounds(r_in,num_particles);
    slice_bunch(r_in,num_particles,nslice,bounds,weight,xpos,ypos,zpos,countslc,pslice);

    for(i=0;i<nslice;i++){        
        register double pos0 = zpos[i];
        if(countslc[i]>0.0){
          for (ii=0;ii<nslice;ii++){
              register double posi = zpos[ii];
              double ds = posi-pos0;
              if(countslc[ii]>0.0 && -ds>waketableT[0] && -ds<waketableT[nelem-1]){
                register double wi = weight[ii];
                register double dx = xpos[ii];
                register double dy = ypos[ii];
                int index = binarySearch(waketableT,-ds,nelem,0,0);              
                double fieldx = getWake(waketableDX,waketableT,-ds,index);
                double fieldy = getWake(waketableDY,waketableT,-ds,index);
                double fieldx2 = getWake(waketableQX,waketableT,-ds,index);
                double fieldy2 = getWake(waketableQY,waketableT,-ds,index);
                double fieldz = getWake(waketableZ,waketableT,-ds,index);
                kx[i] += fx*wi*fieldx*dx;
                ky[i] += fy*wi*fieldy*dy;
                kx2[i] += fqx*wi*fieldx2;
                ky2[i] += fqy*wi*fieldy2;
                kz[i] += fz*wi*fieldz;
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
        long nslice,nelem;
        double on_x,on_y,on_qx,on_qy,on_z;
        double intensity, wakefact, normfactx,normfacty;
        double *waketableT;
        double *waketableDX;
        double *waketableDY;
        double *waketableQX;
        double *waketableQY;
        double *waketableZ;
        
        nslice=atGetLong(ElemData,"Nslice"); check_error();
        nelem=atGetLong(ElemData,"Nelem"); check_error();
        on_x=atGetDouble(ElemData,"On_x"); check_error();
        on_y=atGetDouble(ElemData,"On_y"); check_error();
        on_qx=atGetDouble(ElemData,"On_qx"); check_error();
        on_qy=atGetDouble(ElemData,"On_qy"); check_error();
        on_z=atGetDouble(ElemData,"On_z"); check_error();
        intensity=atGetDouble(ElemData,"Intensity"); check_error();
        wakefact=atGetDouble(ElemData,"Wakefact"); check_error();
        normfactx=atGetDouble(ElemData,"Normfactx"); check_error();
        normfacty=atGetDouble(ElemData,"Normfacty"); check_error();
        waketableT=atGetDoubleArray(ElemData,"WakeT"); check_error();
        waketableDX=atGetDoubleArray(ElemData,"WakeDX"); check_error();
        waketableDY=atGetDoubleArray(ElemData,"WakeDY"); check_error();
        waketableQX=atGetDoubleArray(ElemData,"WakeQX"); check_error();
        waketableQY=atGetDoubleArray(ElemData,"WakeQY"); check_error();
        waketableZ=atGetDoubleArray(ElemData,"WakeZ"); check_error();
        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->nslice=nslice;
        Elem->nelem=nelem;
        Elem->fx=intensity*wakefact*normfactx*on_x;
        Elem->fy=intensity*wakefact*normfacty*on_y;
        Elem->fqx=intensity*wakefact*normfactx*on_qx;
        Elem->fqy=intensity*wakefact*normfacty*on_qy;
        Elem->fz=intensity*wakefact*on_z;
        Elem->waketableT=waketableT;
        Elem->waketableDX=waketableDX;
        Elem->waketableDY=waketableDY;
        Elem->waketableQX=waketableQX;
        Elem->waketableQY=waketableQY;
        Elem->waketableZ=waketableZ;
    }
    impedance_tablePass(r_in,num_particles,Elem);
    return Elem;
}

MODULE_DEF(ImpedanceTablePass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        struct elem El, *Elem=&El;

        long nslice,nelem;
        double on_x,on_y,on_qx,on_qy,on_z;
        double intensity, wakefact, normfactx,normfacty;
        double *waketableT;
        double *waketableDX;
        double *waketableDY;
        double *waketableQX;
        double *waketableQY;
        double *waketableZ;
        
        nslice=atGetLong(ElemData,"Nslice"); check_error();
        nelem=atGetLong(ElemData,"Nelem"); check_error();
        on_x=atGetDouble(ElemData,"On_x"); check_error();
        on_y=atGetDouble(ElemData,"On_y"); check_error();
        on_qx=atGetDouble(ElemData,"On_qx"); check_error();
        on_qy=atGetDouble(ElemData,"On_qy"); check_error();
        on_z=atGetDouble(ElemData,"On_z"); check_error();
        intensity=atGetDouble(ElemData,"Intensity"); check_error();
        wakefact=atGetDouble(ElemData,"Wakefact"); check_error();
        normfactx=atGetDouble(ElemData,"Normfactx"); check_error();
        normfacty=atGetDouble(ElemData,"Normfacty"); check_error();
        waketableT=atGetDoubleArray(ElemData,"WakeT"); check_error();
        waketableDX=atGetDoubleArray(ElemData,"WakeDX"); check_error();
        waketableDY=atGetDoubleArray(ElemData,"WakeDY"); check_error();
        waketableQX=atGetDoubleArray(ElemData,"WakeQX"); check_error();
        waketableQY=atGetDoubleArray(ElemData,"WakeQY"); check_error();
        waketableZ=atGetDoubleArray(ElemData,"WakeZ"); check_error();
        
        Elem->nslice=nslice;
        Elem->nelem=nelem;
        Elem->fx=intensity*wakefact*normfactx*on_x;
        Elem->fy=intensity*wakefact*normfacty*on_y;
        Elem->fqx=intensity*wakefact*normfactx*on_qx;
        Elem->fqy=intensity*wakefact*normfacty*on_qy;
        Elem->fz=intensity*wakefact*on_z;
        Elem->waketableT=waketableT;
        Elem->waketableDX=waketableDX;
        Elem->waketableDY=waketableDY;
        Elem->waketableQX=waketableQX;
        Elem->waketableQY=waketableQY;
        Elem->waketableZ=waketableZ;

        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix: particle array");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        impedance_tablePass(r_in, num_particles, Elem);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(17,1);
        mxSetCell(plhs[0],0,mxCreateString("Nslice"));
        mxSetCell(plhs[0],1,mxCreateString("Intensity"));
        mxSetCell(plhs[0],2,mxCreateString("Wakefact"));
        mxSetCell(plhs[0],3,mxCreateString("WakeT"));
        mxSetCell(plhs[0],4,mxCreateString("WakeDX"));
        mxSetCell(plhs[0],5,mxCreateString("WakeDY"));
        mxSetCell(plhs[0],6,mxCreateString("WakeQX"));
        mxSetCell(plhs[0],7,mxCreateString("WakeQY"));
        mxSetCell(plhs[0],8,mxCreateString("WakeZ"));
        mxSetCell(plhs[0],9,mxCreateString("Nelem"));
        mxSetCell(plhs[0],10,mxCreateString("Normfactx"));
        mxSetCell(plhs[0],11,mxCreateString("Normfacty"));
        mxSetCell(plhs[0],12,mxCreateString("On_x"));
        mxSetCell(plhs[0],13,mxCreateString("On_y"));
        mxSetCell(plhs[0],14,mxCreateString("On_qx"));
        mxSetCell(plhs[0],15,mxCreateString("On_qy"));
        mxSetCell(plhs[0],16,mxCreateString("On_z"));

        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
        }
    }else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif
