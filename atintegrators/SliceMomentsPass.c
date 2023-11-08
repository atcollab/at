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

    int i,ii,ib;
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

    double *xpos = means + turn*nslice*nbunch;
    double *ypos = means + (2*nturns + turn)*nslice*nbunch;
    double *dppos = means + (4*nturns + turn)*nslice*nbunch;
    double *zpos = means + (5*nturns + turn)*nslice*nbunch;
    
    double *xstd = stds + turn*nslice*nbunch;
    double *ystd = stds + (2*nturns + turn)*nslice*nbunch;
    double *dpstd = stds + (4*nturns + turn)*nslice*nbunch;
    double *zstd = stds + (5*nturns + turn)*nslice*nbunch;
    
    double *weight = weights + turn*nslice*nbunch;
    

    /*slices sorted from head to tail (increasing ct)*/
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
            xpos[ii] += rtmp[0];
            ypos[ii] += rtmp[2];
            dppos[ii] += rtmp[4];
            zpos[ii] += rtmp[5];
            xstd[ii] += rtmp[0]*rtmp[0];
            ystd[ii] += rtmp[2]*rtmp[2];
            dpstd[ii] += rtmp[4]*rtmp[4];
            zstd[ii] += rtmp[5]*rtmp[5];
        }
    }

    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE,np_bunch,nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,xpos,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ypos,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,dppos,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,zpos,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    MPI_Allreduce(MPI_IN_PLACE,xstd,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ystd,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,dpstd,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,zstd,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    MPI_Allreduce(MPI_IN_PLACE,weight,nslice*nbunch,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    /*Compute average x/y position and weight of each slice */
    for (i=0;i<nslice*nbunch;i++) {
        ib = (int)(i/nslice);
        if (weight[i] > 0){
            xpos[i] = xpos[i]/weight[i];
            ypos[i] = ypos[i]/weight[i];
            dppos[i] = dppos[i]/weight[i];
            zpos[i] = zpos[i]/weight[i];
            xstd[i] = sqrt(xstd[i]/weight[i]-xpos[i]*xpos[i]);
            ystd[i] = sqrt(ystd[i]/weight[i]-ypos[i]*ypos[i]);
            dpstd[i] = sqrt(dpstd[i]/weight[i]-dppos[i]*dppos[i]);
            zstd[i] = sqrt(zstd[i]/weight[i]-zpos[i]*zpos[i]);    
        }else{
            xpos[i] = 0.0;
            ypos[i] = 0.0;
            dppos[i] = 0.0;
            zpos[i] = 0.0;
            xstd[i] = 0.0;
            ystd[i] = 0.0;
            dpstd[i] = 0.0;
            zstd[i] = smin[ib]+(i%nslice+0.5)*hz[ib];
        }
        zpos[i] += bunch_spos[ib]-bunch_spos[nbunch-1];
        weight[i] *= bunch_currents[ib]/np_bunch[ib];
    }
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
    
    if((turn>=startturn) && (turn<=endturn)){
        slice_beam(r_in, num_particles, nslice, turn-startturn, nturns, nbunch,
                   weights, means, stds, z_cuts, bunch_currents, bunch_spos);
    }; 
};


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                                      double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double *means, *stds, *weights;
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
        int dims[] = {4, Param->nbunch*nslice, endturn-startturn};
        int dimsw[] = {Param->nbunch*nslice,  endturn-startturn};        
        means = atGetDoubleArray(ElemData,"_means"); check_error();
        stds = atGetDoubleArray(ElemData,"_stds"); check_error();
        weights = atGetDoubleArray(ElemData,"_weights"); check_error();
        atCheckArrayDims(ElemData,"_means", 3, dims); check_error();
        atCheckArrayDims(ElemData,"_stds", 3, dims); check_error();
        atCheckArrayDims(ElemData,"_weights", 2, dimsw); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->stds = stds;
        Elem->means = means;
        Elem->weights = weights;
        Elem->turn = 0;
        Elem->startturn = startturn;
        Elem->endturn = endturn;
        Elem->nslice = nslice;
    }   
    SliceMomentsPass(r_in, Param->nbunch, Param->bunch_spos,
                     Param->bunch_currents, num_particles, Elem);
    Elem->turn++;
    return Elem;
}

MODULE_DEF(BeamMomentsPass)        /* Dummy module initialisation */
#endif

#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
    
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        struct elem El, *Elem=&El;
        
        double *means, *stds, weights;
        int startturn = atGetLong(ElemData,"startturn"); check_error();
        int endturn = atGetLong(ElemData,"endturn"); check_error();
        int nslice = atGetLong(ElemData,"nslice"); check_error();    
        means = atGetDoubleArray(ElemData,"_means"); check_error();
        stds = atGetDoubleArray(ElemData,"_stds"); check_error();
        weights = atGetDoubleArray(ElemData,"_weights"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->stds = stds;
        Elem->means = means;
        Elem->weights = weigts;
        Elem->turn = 0;
        Elem->startturn = startturn;
        Elem->endturn = endturn;
        Elem->nslice = nslice;
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
        plhs[0] = mxCreateCellMatrix(2,1);
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
