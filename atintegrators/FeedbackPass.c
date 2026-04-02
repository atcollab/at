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
  double GX;
  double GY;
  double GZ;
  double *closed_orbit;
  double *xbuffer;
  double *ybuffer;
  double *zbuffer;
  long bufferlength_x;
  long bufferlength_y;
  long bufferlength_z;
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
    double mx[] = {0.0};
    double mx_set[1] = {0.0};
    double my[] = {0.0};
    double my_set[1] = {0.0};
    double mz[] = {0.0};
    double mz_set[1] = {0.0};
    long npart = 0.0;
            
    double *closed_orbit; 
    closed_orbit = Elem->closed_orbit;
    
    double gx = Elem->GX;
    double gy = Elem->GY;
    double gz = Elem->GZ;
    
    double *xbuffer = Elem->xbuffer;
    double *ybuffer = Elem->ybuffer;
    double *zbuffer = Elem->zbuffer;
    long bufferlength_x = Elem->bufferlength_x;
    long bufferlength_y = Elem->bufferlength_y;
    long bufferlength_z = Elem->bufferlength_z;
    
    long bufferlengthnow_x = 0;
    long bufferlengthnow_y = 0;
    long bufferlengthnow_z = 0;
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        double *r6 = r_in+c*6;

        if (!atIsNaN(r6[0])) {
            npart += 1; 
            mx[0] += r6[0];
            my[0] += r6[2];
            mz[0] += r6[5]; 
        }
    }

    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE, &npart, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &mx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(MPI_IN_PLACE, &my, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(MPI_IN_PLACE, &mz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
    
    if (npart>0.0){
        mx[0] /= npart;
        my[0] /= npart;
        mz[0] /= npart;
    }

    // horizonal
    if(bufferlength_x>0){
        write_buffer(mx, xbuffer, 1, bufferlength_x);
        bufferlengthnow_x = check_buffer_length(xbuffer, bufferlength_x, 1);
        if(bufferlengthnow_x == bufferlength_x){
            compute_buffer_mean(mx_set, xbuffer, bufferlength_x, bufferlength_x, 1);
        }
        
    }else{
        mx_set[0] = mx[0];
    }

    // vertical    
    if(bufferlength_y>0){
        write_buffer(my, ybuffer, 1, bufferlength_y);
        bufferlengthnow_y = check_buffer_length(ybuffer, bufferlength_y, 1);
        if(bufferlengthnow_y == bufferlength_y){
           compute_buffer_mean(my_set, ybuffer, bufferlength_y, bufferlength_y, 1);
        }
        
    }else{
        my_set[0] = my[0];
    }
    
    // longitudinal
    if(bufferlength_z>0){
        write_buffer(mz, zbuffer, 1, bufferlength_z);
        bufferlengthnow_z = check_buffer_length(zbuffer, bufferlength_z, 1);
        if(bufferlengthnow_z == bufferlength_z){
            compute_buffer_mean(mz_set, zbuffer, bufferlength_z, bufferlength_z, 1);
        }

    }else{
        mz_set[0] = mz[0];
    }
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        double *r6 = r_in+c*6;
        if (!atIsNaN(r6[0])) {
            r6[0] -= (mx_set[0]-closed_orbit[0])*gx;
            r6[2] -= (my_set[0]-closed_orbit[2])*gy;
            r6[5] -= (mz_set[0]-closed_orbit[5])*gz;    
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double GX, GY, GZ, *closed_orbit;
        double *xbuffer, *ybuffer, *zbuffer;
        long bufferlength_x, bufferlength_y, bufferlength_z;
        GX=atGetOptionalDouble(ElemData,"Gx",0.0); check_error();
        GY=atGetOptionalDouble(ElemData,"Gy",0.0); check_error();
        GZ=atGetOptionalDouble(ElemData,"Gz",0.0); check_error();
        closed_orbit = atGetDoubleArray(ElemData,"closed_orbit"); check_error();
        xbuffer=atGetDoubleArray(ElemData,"_xbuffer"); check_error();
        ybuffer=atGetDoubleArray(ElemData,"_ybuffer"); check_error();
        zbuffer=atGetDoubleArray(ElemData,"_zbuffer"); check_error();
        bufferlength_x=atGetLong(ElemData,"_bufferlength_x"); check_error();
        bufferlength_y=atGetLong(ElemData,"_bufferlength_y"); check_error();        
        bufferlength_z=atGetLong(ElemData,"_bufferlength_z"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->GX=GX;
        Elem->GY=GY;
        Elem->GZ=GZ;
        Elem->closed_orbit = closed_orbit;
        Elem->xbuffer = xbuffer;
        Elem->ybuffer = ybuffer;
        Elem->zbuffer = zbuffer;
        Elem->bufferlength_x = bufferlength_x;
        Elem->bufferlength_y = bufferlength_y;
        Elem->bufferlength_z = bufferlength_z;

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
        
        double GX,GY,GZ,*closed_orbit;
        double *xbuffer, *ybuffer, *zbuffer;
        long bufferlength_x, bufferlength_y, bufferlength_z;
        GX=atGetOptionalDouble(ElemData,"Gx",0.0); check_error();
        GY=atGetOptionalDouble(ElemData,"Gy",0.0); check_error();
        GZ=atGetOptionalDouble(ElemData,"Gz",0.0); check_error();
        closed_orbit = atGetDoubleArray(ElemData,"closed_orbit");check_error();
        xbuffer=atGetDoubleArray(ElemData,"_xbuffer"); check_error();
        ybuffer=atGetDoubleArray(ElemData,"_ybuffer"); check_error();
        zbuffer=atGetDoubleArray(ElemData,"_zbuffer"); check_error();
        bufferlength_x=atGetLong(ElemData,"_bufferlength_x"); check_error();
        bufferlength_y=atGetLong(ElemData,"_bufferlength_y"); check_error();
        bufferlength_z=atGetLong(ElemData,"_bufferlength_z"); check_error();
        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->GX=GX;
        Elem->GY=GY;
        Elem->GZ=GZ;
        Elem->closed_orbit = closed_orbit;
        Elem->xbuffer = xbuffer;
        Elem->ybuffer = ybuffer;
        Elem->zbuffer = zbuffer;
        Elem->bufferlength_x = bufferlength_x;
        Elem->bufferlength_y = bufferlength_y;
        Elem->bufferlength_z = bufferlength_z;
        
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
        mxSetCell(plhs[0],1,mxCreateString("_xbuffer"));
        mxSetCell(plhs[0],2,mxCreateString("_ybuffer"));        
        mxSetCell(plhs[0],3,mxCreateString("_zbuffer"));        
        mxSetCell(plhs[0],4,mxCreateString("_bufferlength_x"));                
        mxSetCell(plhs[0],5,mxCreateString("_bufferlength_y"));                
        mxSetCell(plhs[0],6,mxCreateString("_bufferlength_z"));                

        if(nlhs>1) /* optional fields */
        {
            plhs[1] = mxCreateCellMatrix(3,1);
            mxSetCell(plhs[1],0,mxCreateString("Gx"));
            mxSetCell(plhs[1],1,mxCreateString("Gy"));
            mxSetCell(plhs[1],2,mxCreateString("Gz"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif

