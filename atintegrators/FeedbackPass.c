/* IdentityPass.c 
   Accelerator Toolbox
   Revision 7/16/03
   A.Terebilo
*/

#include "atelem.c"
#include "atlalib.c"
#ifdef MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif

struct elem 
{
  double Gxp;
  double Gyp;
  double Gdp;
  double *closed_orbit;
  double *xpbuffer;
  double *ypbuffer;
  double *dpbuffer;
  long bufferlength_xp;
  long bufferlength_yp;
  long bufferlength_dp;
};

void write_buffer(double *data, double *buffer, int datasize, int buffersize){
    if(buffersize>1){
        memmove(buffer, buffer + datasize, datasize*(buffersize-1)*sizeof(double));
    }
    memcpy(buffer + datasize*(buffersize-1), data, datasize*sizeof(double));
}

int check_buffer_length(double *buffer, long buffersize, long numcolumns){
    int c;
    int bufferlengthnow=0;
    for (c=0; c<numcolumns*buffersize; c++){
        if (buffer[c]!=0.0){
            bufferlengthnow += 1;
        }
    }
    bufferlengthnow /= numcolumns;
    return bufferlengthnow;
}

static void compute_buffer_mean(double *out_array, double *buffer, long windowlength, long buffersize, long numcolumns){

    int c,p,offset;
    offset = buffersize - windowlength;

    for (p=0; p<numcolumns; p++) {
        out_array[p] = 0.0;
    }
    
    for (c=offset; c<buffersize; c++) {
        for (p=0; p<numcolumns; p++) {
            out_array[p] += buffer[numcolumns*c+p];
        }
    }
    
    for (p=0; p<numcolumns; p++) {
        out_array[p] /= windowlength ; 
    }
}



void FeedbackPass(double *r_in, int num_particles, struct elem *Elem)
{	
    int c;
    double mxp[] = {0.0}; // the mean that will be computed
    double mxp_set[1] = {0.0}; // define pointer for the mean that will be set (buffer or one turn)
    double myp[] = {0.0};
    double myp_set[1] = {0.0};
    double mdp[] = {0.0};
    double mdp_set[1] = {0.0};
    long npart = 0.0;
            
    double *closed_orbit; 
    closed_orbit = Elem->closed_orbit;
    
    double gxp = Elem->Gxp;
    double gyp = Elem->Gyp;
    double gdp = Elem->Gdp;
    
    double *xpbuffer = Elem->xpbuffer;
    double *ypbuffer = Elem->ypbuffer;
    double *dpbuffer = Elem->dpbuffer;
    long bufferlength_xp = Elem->bufferlength_xp;
    long bufferlength_yp = Elem->bufferlength_yp;
    long bufferlength_dp = Elem->bufferlength_dp;
    
    long bufferlengthnow_xp = 0;
    long bufferlengthnow_yp = 0;
    long bufferlengthnow_dp = 0;
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        double *r6 = r_in+c*6;

        if (!atIsNaN(r6[0])) {
            npart += 1; 
            mxp[0] += r6[1]; //xp
            myp[0] += r6[3]; //yp
            mdp[0] += r6[4]; //dp
        }
    }

    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE, &npart, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &mxp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(MPI_IN_PLACE, &myp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(MPI_IN_PLACE, &mdp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
    
    if (npart>0.0){
        mxp[0] /= npart;
        myp[0] /= npart;
        mdp[0] /= npart;
    }

    // horizonal
    if(bufferlength_xp>0){
        write_buffer(mxp, xpbuffer, 1, bufferlength_xp);
        bufferlengthnow_xp = check_buffer_length(xpbuffer, bufferlength_xp, 1);
        if(bufferlengthnow_xp == bufferlength_xp){
            compute_buffer_mean(mxp_set, xpbuffer, bufferlength_xp, bufferlength_xp, 1);
        }
        
    }else{
        mxp_set[0] = mxp[0];
    }

    // vertical    
    if(bufferlength_yp>0){
        write_buffer(myp, ypbuffer, 1, bufferlength_yp);
        bufferlengthnow_yp = check_buffer_length(ypbuffer, bufferlength_yp, 1);
        if(bufferlengthnow_yp == bufferlength_yp){
           compute_buffer_mean(myp_set, ypbuffer, bufferlength_yp, bufferlength_yp, 1);
        }
        
    }else{
        myp_set[0] = myp[0];
    }
    
    // longitudinal
    if(bufferlength_dp>0){
        write_buffer(mdp, dpbuffer, 1, bufferlength_dp);
        bufferlengthnow_dp = check_buffer_length(dpbuffer, bufferlength_dp, 1);
        if(bufferlengthnow_dp == bufferlength_dp){
            compute_buffer_mean(mdp_set, dpbuffer, bufferlength_dp, bufferlength_dp, 1);
        }

    }else{
        mdp_set[0] = mdp[0];
    }
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        double *r6 = r_in+c*6;
        if (!atIsNaN(r6[0])) {
            r6[1] -= (mxp_set[0]-closed_orbit[0])*gxp;
            r6[3] -= (myp_set[0]-closed_orbit[1])*gyp;
            r6[4] -= (mdp_set[0]-closed_orbit[2])*gdp;    
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Gxp, Gyp, Gdp, *closed_orbit;
        double *xpbuffer, *ypbuffer, *dpbuffer;
        long bufferlength_xp, bufferlength_yp, bufferlength_dp;
        Gxp=atGetOptionalDouble(ElemData,"Gxp",0.0); check_error();
        Gyp=atGetOptionalDouble(ElemData,"Gyp",0.0); check_error();
        Gdp=atGetOptionalDouble(ElemData,"Gdp",0.0); check_error();
        closed_orbit = atGetDoubleArray(ElemData,"closed_orbit"); check_error();
        xpbuffer=atGetDoubleArray(ElemData,"_xpbuffer"); check_error();
        ypbuffer=atGetDoubleArray(ElemData,"_ypbuffer"); check_error();
        dpbuffer=atGetDoubleArray(ElemData,"_dpbuffer"); check_error();
        bufferlength_xp=atGetLong(ElemData,"_bufferlength_xp"); check_error();
        bufferlength_yp=atGetLong(ElemData,"_bufferlength_yp"); check_error();        
        bufferlength_dp=atGetLong(ElemData,"_bufferlength_dp"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Gxp=Gxp;
        Elem->Gyp=Gyp;
        Elem->Gdp=Gdp;
        Elem->closed_orbit = closed_orbit;
        Elem->xpbuffer = xpbuffer;
        Elem->ypbuffer = ypbuffer;
        Elem->dpbuffer = dpbuffer;
        Elem->bufferlength_xp = bufferlength_xp;
        Elem->bufferlength_yp = bufferlength_yp;
        Elem->bufferlength_dp = bufferlength_dp;

    }
    FeedbackPass(r_in,num_particles,Elem);
    return Elem;
}

MODULE_DEF(FeedbackPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/


#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
    
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        struct elem El, *Elem=&El;
        
        double Gxp,Gyp,Gdp,*closed_orbit;
        double *xpbuffer, *ypbuffer, *dpbuffer;
        long bufferlength_xp, bufferlength_yp, bufferlength_dp;
        Gxp=atGetOptionalDouble(ElemData,"Gxp",0.0); check_error();
        Gyp=atGetOptionalDouble(ElemData,"Gyp",0.0); check_error();
        Gdp=atGetOptionalDouble(ElemData,"Gdp",0.0); check_error();
        closed_orbit = atGetDoubleArray(ElemData,"closed_orbit");check_error();
        xpbuffer=atGetDoubleArray(ElemData,"_xpbuffer"); check_error();
        ypbuffer=atGetDoubleArray(ElemData,"_ypbuffer"); check_error();
        dpbuffer=atGetDoubleArray(ElemData,"_dpbuffer"); check_error();
        bufferlength_xp=atGetLong(ElemData,"_bufferlength_xp"); check_error();
        bufferlength_yp=atGetLong(ElemData,"_bufferlength_yp"); check_error();
        bufferlength_dp=atGetLong(ElemData,"_bufferlength_dp"); check_error();
        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Gxp=Gxp;
        Elem->Gyp=Gyp;
        Elem->Gdp=Gdp;
        Elem->closed_orbit = closed_orbit;
        Elem->xpbuffer = xpbuffer;
        Elem->ypbuffer = ypbuffer;
        Elem->dpbuffer = dpbuffer;
        Elem->bufferlength_xp = bufferlength_xp;
        Elem->bufferlength_yp = bufferlength_yp;
        Elem->bufferlength_dp = bufferlength_dp;
        
    if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix: particle array");
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        FeedbackPass(r_in,num_particles,Elem);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(7,1);
        mxSetCell(plhs[0],0,mxCreateString("closed_orbit"));
        mxSetCell(plhs[0],1,mxCreateString("_xpbuffer"));
        mxSetCell(plhs[0],2,mxCreateString("_ypbuffer"));        
        mxSetCell(plhs[0],3,mxCreateString("_dpbuffer"));        
        mxSetCell(plhs[0],4,mxCreateString("_bufferlength_xp"));                
        mxSetCell(plhs[0],5,mxCreateString("_bufferlength_yp"));                
        mxSetCell(plhs[0],6,mxCreateString("_bufferlength_dp"));                

        if(nlhs>1) /* optional fields */
        {
            plhs[1] = mxCreateCellMatrix(3,1);
            mxSetCell(plhs[1],0,mxCreateString("Gxp"));
            mxSetCell(plhs[1],1,mxCreateString("Gyp"));
            mxSetCell(plhs[1],2,mxCreateString("Gdp"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif

