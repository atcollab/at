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
  double *sposs;
  double *weights;
  double *z_cuts;
};


static void slice_beam(double *r_in,int num_particles,int nslice,int turn,
                       int nturns, int nbunch, double *weights, double *sposs,
                       double *means, double *stds, double *z_cuts,
                       double *bunch_currents, double beam_current){

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
    
    void *buffer = atCalloc(8*nbunch*nslice, sizeof(double));
    double *dptr = (double *) buffer;
    int idx[] = {0, 2, 4};
    double *pos = dptr; dptr += 3*nbunch*nslice;
    double *std = dptr; dptr += 3*nbunch*nslice;
    double *spos = dptr; dptr += nbunch*nslice;
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
            for(iii=0; iii<3; iii++) {
                pos[iii+ii*3] += rtmp[idx[iii]];
                std[iii+ii*3] += rtmp[idx[iii]]*rtmp[idx[iii]];
            }
        }
    }

    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE,np_bunch,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,pos,3*nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,std,3*nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,weight,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    for (i=0;i<nslice*nbunch;i++) {
        ib = (int)(i/nslice);
        for(ii=0; ii<3; ii++){
            if (weight[i] > 0){
                pos[3*i+ii] = pos[3*i+ii]/weight[i];
                std[3*i+ii] = sqrt(std[3*i+ii]/weight[i]-pos[3*i+ii]*pos[3*i+ii]);
            }
            else{
                pos[3*i+ii] = NAN;
                std[3*i+ii] = NAN;
            }
        }
        spos[i] = smin[ib]+(i%nslice+0.5)*hz[ib];
        if(beam_current>0.0){
            weight[i] *= bunch_currents[ib]/np_bunch[ib];
        }else{
            weight[i] *= 1.0/np_bunch[ib]; 
        }    
    }
    
    means += 3*nbunch*nslice*turn;
    stds += 3*nbunch*nslice*turn;
    sposs += nbunch*nslice*turn;
    weights += nbunch*nslice*turn;
    memcpy(means, pos, 3*nbunch*nslice*sizeof(double));
    memcpy(stds, std, 3*nbunch*nslice*sizeof(double));
    memcpy(sposs, spos, nbunch*nslice*sizeof(double));
    memcpy(weights, weight, nbunch*nslice*sizeof(double));
    
    atFree(buffer);
    atFree(np_bunch);
    atFree(smin);
    atFree(smax);
    atFree(hz);
};


void SliceMomentsPass(double *r_in, int nbunch, double *bunch_currents,
                      double beam_current, int num_particles, struct elem *Elem) {

    int startturn = Elem->startturn;
    int endturn = Elem->endturn;
    int nturns = endturn-startturn;
    int turn = Elem->turn;
    int nslice = Elem->nslice;
    double *stds = Elem->stds;
    double *means = Elem->means;
    double *sposs = Elem->sposs;
    double *weights = Elem->weights;
    double *z_cuts = Elem->z_cuts;
    
    if((turn>=startturn) && (turn<endturn)){
        slice_beam(r_in, num_particles, nslice, turn-startturn, nturns, nbunch,
                   weights, sposs, means, stds, z_cuts, bunch_currents, beam_current);
    }; 
};


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                                      double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double *means, *stds, *sposs, *weights, *z_cuts;
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
        int dims[] = {3, Param->nbunch*nslice, endturn-startturn};
        int dimsw[] = {Param->nbunch*nslice,  endturn-startturn};
        means = atGetDoubleArray(ElemData,"_means"); check_error();
        stds = atGetDoubleArray(ElemData,"_stds"); check_error();
        sposs = atGetDoubleArray(ElemData,"_spos"); check_error();
        weights = atGetDoubleArray(ElemData,"_weights"); check_error();
        z_cuts=atGetOptionalDoubleArray(ElemData,"ZCuts"); check_error();
        atCheckArrayDims(ElemData,"_means", 3, dims); check_error();
        atCheckArrayDims(ElemData,"_stds", 3, dims); check_error();
        atCheckArrayDims(ElemData,"_spos", 2, dimsw); check_error();
        atCheckArrayDims(ElemData,"_weights", 2, dimsw); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->stds = stds;
        Elem->means = means;
        Elem->sposs = sposs;
        Elem->weights = weights;
        Elem->turn = 0;
        Elem->startturn = startturn;
        Elem->endturn = endturn;
        Elem->nslice = nslice;
        Elem->z_cuts = z_cuts;
    }   
    SliceMomentsPass(r_in, Param->nbunch, Param->bunch_currents,
                     Param->beam_current, num_particles, Elem);
    Elem->turn++;
    return Elem;
}

MODULE_DEF(SliceMomentsPass)        /* Dummy module initialisation */
#endif

#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
    
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
        double *bcurr = malloc(sizeof(double));
        bcurr[0] = 0.0;
        SliceMomentsPass(r_in,1,bcurr, 1.0,num_particles,Elem);
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
