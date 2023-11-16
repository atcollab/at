#include "atelem.c"
#include "atimplib.c"
#include <math.h>
#include <float.h>
#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif


struct elem
{
  int startturn;
  int endturn;
  int turn;
  int nslice;
  double *stds;
  double *means;
  double *weights;
  double *z_cuts;
};


static void slice_beam(double *r_in,int num_particles,int nslice,int turn,
                       int nturns, int nbunch, double *weights, double *means,
                       double *stds, double *z_cuts, double *bunch_currents,
                       double *bunch_spos){

    int i,ii,iii,ib;
    double *rtmp;

    double *smin = atMalloc(nbunch*sizeof(double));
    double *smax = atMalloc(nbunch*sizeof(double));
    double *hz = atMalloc(nbunch*sizeof(double));
    double *np_bunch = atMalloc(nbunch*sizeof(double));
    getbounds(r_in,nbunch,num_particles,smin,smax,z_cuts);

    for(i=0;i<nbunch;i++){
        hz[i] = (smax[i]-smin[i])/(nslice);
        np_bunch[i] = 0.0;
    }
    
    void *buffer = atCalloc(9*nbunch*nslice, sizeof(double));
    double *dptr = (double *) buffer;
    int idx[] = {0, 2, 4, 5};    
    double *pos = dptr; dptr += 4*nbunch*nslice;
    double *std = dptr; dptr += 4*nbunch*nslice;
    double *weight = dptr;
    
    for (i=0;i<num_particles;i++) {
        rtmp = r_in+i*6;
        ib = i%nbunch;
        np_bunch[ib] += 1.0;
        if (!atIsNaN(rtmp[0]) && (rtmp[5] >= smin[ib]) && (rtmp[5] <= smax[ib])) {
            if (rtmp[5] == smax[ib]){
                ii = nslice-1 + ib*nslice;
            }
            else {
                ii = (int)(floor((rtmp[5]-smin[ib])/hz[ib])) + ib*nslice;
            }
            weight[ii] += 1.0;
            for(iii=0; iii<4; iii++) {
                pos[iii+ii*4] += rtmp[idx[iii]];
                std[iii+ii*4] += rtmp[idx[iii]]*rtmp[idx[iii]];
            }
        }
    }

    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE,np_bunch,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,pos,4*nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);   
    MPI_Allreduce(MPI_IN_PLACE,std,4*nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);   
    MPI_Allreduce(MPI_IN_PLACE,weight,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    for (i=0;i<nslice*nbunch;i++) {
        ib = (int)(i/nslice);
        for(ii=0; ii<4; ii++){
            if (weight[i] > 0){
                pos[4*i+ii] = pos[4*i+ii]/weight[i];
                std[4*i+ii] = sqrt(std[4*i+ii]/weight[i]-pos[4*i+ii]*pos[4*i+ii]);
            }   
            else{
                pos[4*i+ii] = (ii==3) ? smin[ib]+(i%nslice+0.5)*hz[ib] : NAN;
                std[4*i+ii] = NAN;
            }                 
        }
        pos[i+4] += bunch_spos[ib]-bunch_spos[nbunch-1];
        weight[i] *= bunch_currents[ib]/np_bunch[ib];
    }
    
    means += 4*nbunch*nslice*turn;
    stds += 4*nbunch*nslice*turn;
    weights += nbunch*nslice*turn;
    memcpy(means, pos, 4*nbunch*nslice*sizeof(double));
    memcpy(stds, std, 4*nbunch*nslice*sizeof(double));
    memcpy(weights, weight, nbunch*nslice*sizeof(double));
    
    atFree(buffer);
    atFree(np_bunch);
    atFree(smin);
    atFree(smax);
    atFree(hz);
};


void SliceMomentsPass(double *r_in, int nbunch, double *bunch_spos, double *bunch_currents,
                      int num_particles, struct elem *Elem) {

    int startturn = Elem->startturn;
    int endturn = Elem->endturn;
    int nturns = endturn-startturn;
    int turn = Elem->turn;
    int nslice = Elem->nslice;
    double *stds = Elem->stds;
    double *means = Elem->means;
    double *weights = Elem->weights;
    double *z_cuts = Elem->z_cuts;
    
    if((turn>=startturn) && (turn<endturn)){
        slice_beam(r_in, num_particles, nslice, turn-startturn, nturns, nbunch,
                   weights, means, stds, z_cuts, bunch_currents, bunch_spos);
    }; 
};


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                                      double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double *means, *stds, *weights, *z_cuts;
        int nslice = atGetLong(ElemData,"nslice"); check_error();
        int startturn = atGetLong(ElemData,"_startturn"); check_error();
        int endturn = atGetLong(ElemData,"_endturn"); check_error();
        if (endturn<0 || startturn<0){
            atError("starturn and endturn have to be greater than 0");
        } else if (endturn<0 || startturn<0){
            atError("starturn has to be smaller than endturn.");
        } else if (endturn > Param->num_turns){
            atWarning("endturn exceed the total number of turns");
        };
        int dims[] = {4, Param->nbunch, nslice, endturn-startturn};
        int dimsw[] = {Param->nbunch, nslice,  endturn-startturn};
        means = atGetDoubleArray(ElemData,"_means"); check_error();
        stds = atGetDoubleArray(ElemData,"_stds"); check_error();
        weights = atGetDoubleArray(ElemData,"_weights"); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
        atCheckArrayDims(ElemData,"_means", 4, dims); check_error();
        atCheckArrayDims(ElemData,"_stds", 4, dims); check_error();
        atCheckArrayDims(ElemData,"_weights", 3, dimsw); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->stds = stds;
        Elem->means = means;
        Elem->weights = weights;
        Elem->turn = 0;
        Elem->startturn = startturn;
        Elem->endturn = endturn;
        Elem->nslice = nslice;
        Elem->z_cuts = z_cuts;
    }   
    SliceMomentsPass(r_in, Param->nbunch, Param->bunch_spos,
                     Param->bunch_currents, num_particles, Elem);
    Elem->turn++;
    return Elem;
}

MODULE_DEF(SliceMomentsPass)        /* Dummy module initialisation */
#endif

#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
    
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        struct elem El, *Elem=&El;
        
        double *means, *stds, *weights, *z_cuts;
        int startturn = atGetLong(ElemData,"startturn"); check_error();
        int endturn = atGetLong(ElemData,"endturn"); check_error();
        int nslice = atGetLong(ElemData,"nslice"); check_error();    
        means = atGetDoubleArray(ElemData,"_means"); check_error();
        stds = atGetDoubleArray(ElemData,"_stds"); check_error();
        weights = atGetDoubleArray(ElemData,"_weights"); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->stds = stds;
        Elem->means = means;
        Elem->weights = weights;
        Elem->turn = 0;
        Elem->startturn = startturn;
        Elem->endturn = endturn;
        Elem->nslice = nslice;
        Elem->z_cuts = z_cuts;
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix: particle array");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
                double *bspos = malloc(sizeof(double));
        double *bcurr = malloc(sizeof(double));
        bspos[0] = 0.0;
        bcurr[0] = 0.0;
        SliceMomentsPass(r_in,1,bspos,bcurr,num_particles,Elem);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("_means"));
        mxSetCell(plhs[0],1,mxCreateString("_stds"));
        mxSetCell(plhs[0],2,mxCreateString("_weights"));
        mxSetCell(plhs[0],3,mxCreateString("startturn"));
        mxSetCell(plhs[0],4,mxCreateString("endturn"));
        mxSetCell(plhs[0],5,mxCreateString("nslice"));
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif
