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


void WakeFieldPass(double *r_in,int num_particles,double circumference,int nbunch,
                   double *bunch_spos,double *bunch_currents,struct elem *Elem) {
    /*
     * r_in - 6-by-N matrix of initial conditions reshaped into
     * 1-d array of 6*N elements
     */   
    long nslice = Elem->nslice;
    long nelem = Elem->nelem;
    long nturns = Elem->nturns;
    double *normfact = Elem->normfact;
    double *waketableT = Elem->waketableT;
    double *waketableDX = Elem->waketableDX;
    double *waketableDY = Elem->waketableDY;
    double *waketableQX = Elem->waketableQX;
    double *waketableQY = Elem->waketableQY;
    double *waketableZ = Elem->waketableZ;
    double *turnhistory = Elem->turnhistory;
    double *z_cuts = Elem->z_cuts;    

    size_t sz = 5*nslice*nbunch*sizeof(double) + num_particles*sizeof(int);
    int c;

    int *pslice;
    double *kx;
    double *ky;
    double *kx2;
    double *ky2;
    double *kz;

    void *buffer = atMalloc(sz);
    double *dptr = (double *) buffer;
    int *iptr;

    kx = dptr; dptr += nslice*nbunch;
    ky = dptr; dptr += nslice*nbunch;
    kx2 = dptr; dptr += nslice*nbunch;
    ky2 = dptr; dptr += nslice*nbunch;
    kz = dptr; dptr += nslice*nbunch;

    iptr = (int *) dptr;
    pslice = iptr; iptr += num_particles;

    /*slices beam and compute kick*/
    rotate_table_history(nturns,nslice*nbunch,turnhistory,circumference);
    slice_bunch(r_in,num_particles,nslice,nturns,nbunch,bunch_spos,bunch_currents,
                turnhistory,pslice,z_cuts);
    compute_kicks(nslice*nbunch,nturns,nelem,turnhistory,waketableT,waketableDX,
                  waketableDY,waketableQX,waketableQY,waketableZ,
                  normfact,kx,ky,kx2,ky2,kz);
    
    /*apply kicks*/
    /* OpenMP not efficient. Too much shared data ?
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
    shared(r_in,num_particles,pslice,kx,kx2,ky,ky2,kz) private(c)
    */
    for (c=0; c<num_particles; c++) {
        double *r6 = r_in+c*6;
        int islice=pslice[c];
        if (!atIsNaN(r6[0])) {
            r6[4] += kz[islice];
            r6[1] += (kx[islice]+r6[0]*kx2[islice])*(1+r6[4]);
            r6[3] += (ky[islice]+r6[2]*ky2[islice])*(1+r6[4]);
        }
    }
    atFree(buffer);
}


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        long nslice,nelem,nturns;
        double wakefact;
        static double lnf[3];
        double *normfact;
        double *waketableT;
        double *waketableDX;
        double *waketableDY;
        double *waketableQX;
        double *waketableQY;
        double *waketableZ;
        double *turnhistory;
        double *z_cuts;
        int i;

        nslice=atGetLong(ElemData,"_nslice"); check_error();
        nelem=atGetLong(ElemData,"_nelem"); check_error();
        nturns=atGetLong(ElemData,"_nturns"); check_error();
        wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
        waketableT=atGetDoubleArray(ElemData,"_wakeT"); check_error();
        turnhistory=atGetDoubleArray(ElemData,"_turnhistory"); check_error();
        normfact=atGetDoubleArray(ElemData,"NormFact"); check_error();
        /*optional attributes*/
        waketableDX=atGetOptionalDoubleArray(ElemData,"_wakeDX"); check_error();
        waketableDY=atGetOptionalDoubleArray(ElemData,"_wakeDY"); check_error();
        waketableQX=atGetOptionalDoubleArray(ElemData,"_wakeQX"); check_error();
        waketableQY=atGetOptionalDoubleArray(ElemData,"_wakeQY"); check_error();
        waketableZ=atGetOptionalDoubleArray(ElemData,"_wakeZ"); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();

        int dimsth[] = {Param->nbunch*nslice*nturns, 4};
        atCheckArrayDims(ElemData,"_turnhistory", 2, dimsth); check_error();
        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->nslice=nslice;
        Elem->nelem=nelem;
        Elem->nturns=nturns;
        for(i=0;i<3;i++){
           lnf[i]=normfact[i]*wakefact;
        }
        Elem->normfact=lnf;
        Elem->waketableT=waketableT;
        Elem->waketableDX=waketableDX;
        Elem->waketableDY=waketableDY;
        Elem->waketableQX=waketableQX;
        Elem->waketableQY=waketableQY;
        Elem->waketableZ=waketableZ;
        Elem->turnhistory=turnhistory;
        Elem->z_cuts=z_cuts;
    }
    if(num_particles<Param->nbunch){
        atError("Number of particles has to be greater or equal to the number of bunches.");
    }else if (num_particles%Param->nbunch!=0){
        atWarning("Number of particles not a multiple of the number of bunches: uneven bunch load.");
    }
    WakeFieldPass(r_in,num_particles,Param->RingLength,Param->nbunch,Param->bunch_spos,
                  Param->bunch_currents,Elem);
    return Elem;
}

MODULE_DEF(WakeFieldPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/


#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        int i;
        struct elem El, *Elem=&El;

        long nslice,nelem,nturns;
        double wakefact;
        static double lnf[3];
        double *normfact;
        double *waketableT;
        double *waketableDX;
        double *waketableDY;
        double *waketableQX;
        double *waketableQY;
        double *waketableZ;
        double *turnhistory;
        double *z_cuts;

        nslice=atGetLong(ElemData,"_nslice"); check_error();
        nelem=atGetLong(ElemData,"_nelem"); check_error();
        nturns=atGetLong(ElemData,"_nturns"); check_error();
        wakefact=atGetDouble(ElemData,"_wakefact"); check_error();
        waketableT=atGetDoubleArray(ElemData,"_wakeT"); check_error();
        turnhistory=atGetDoubleArray(ElemData,"_turnhistory"); check_error();
        normfact=atGetDoubleArray(ElemData,"NormFact"); check_error();
        /*optional attributes*/
        waketableDX=atGetOptionalDoubleArray(ElemData,"_wakeDX"); check_error();
        waketableDY=atGetOptionalDoubleArray(ElemData,"_wakeDY"); check_error();
        waketableQX=atGetOptionalDoubleArray(ElemData,"_wakeQX"); check_error();
        waketableQY=atGetOptionalDoubleArray(ElemData,"_wakeQY"); check_error();
        waketableZ=atGetOptionalDoubleArray(ElemData,"_wakeZ"); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
        
        Elem->nslice=nslice;
        Elem->nelem=nelem;
        Elem->nturns=nturns;
        for(i=0;i<3;i++){
           lnf[i]=normfact[i]*wakefact;
        }
        Elem->normfact=lnf;
        Elem->waketableT=waketableT;
        Elem->waketableDX=waketableDX;
        Elem->waketableDY=waketableDY;
        Elem->waketableQX=waketableQX;
        Elem->waketableQY=waketableQY;
        Elem->waketableZ=waketableZ;
        Elem->turnhistory=turnhistory;
        Elem->z_cuts=z_cuts;

        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix: particle array");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        double *bspos = malloc(sizeof(double));
        double *bcurr = malloc(sizeof(double));
        bspos[0] = 0.0;
        bcurr[0] = 0.0;
        WakeFieldPass(r_in,num_particles, 1, 1, bspos, bcurr, Elem);
        free(bspos);
        free(bcurr);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(7,1);
        mxSetCell(plhs[0],0,mxCreateString("_nelem"));
        mxSetCell(plhs[0],1,mxCreateString("_nslice"));
        mxSetCell(plhs[0],2,mxCreateString("_nturns"));
        mxSetCell(plhs[0],3,mxCreateString("_wakefact"));
        mxSetCell(plhs[0],4,mxCreateString("_wakeT"));
        mxSetCell(plhs[0],5,mxCreateString("_turnhistory"));
        mxSetCell(plhs[0],6,mxCreateString("Normfact"));

        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(6,1); /* No optional fields */
            mxSetCell(plhs[0],0,mxCreateString("_wakeDX"));
            mxSetCell(plhs[0],1,mxCreateString("_wakeDY"));
            mxSetCell(plhs[0],2,mxCreateString("_wakeQX"));
            mxSetCell(plhs[0],3,mxCreateString("_wakeQY"));
            mxSetCell(plhs[0],4,mxCreateString("_wakeZ"));
            mxSetCell(plhs[0],5,mxCreateString("ZCuts"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif
