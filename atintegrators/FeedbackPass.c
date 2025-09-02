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
};

void FeedbackPass(double *r_in, int num_particles, struct elem *Elem)
{	
    int c;
    double *mx = atMalloc(sizeof(double));
    double *my = atMalloc(sizeof(double));
    double *mz = atMalloc(sizeof(double));
    long *npart = atMalloc(sizeof(long));
    
    *mx = 0.0; 
    *my = 0.0;
    *mz = 0.0; 
    *npart = 0;
    
    double *closed_orbit; 
    closed_orbit = Elem->closed_orbit;
    
    double gx = Elem->GX;
    double gy = Elem->GY;
    double gz = Elem->GZ;
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        double *r6 = r_in+c*6;

        if (!atIsNaN(r6[0])) {
            *npart += 1; 
            *mx += r6[0];
            *my += r6[2];
            *mz += r6[5]; 
        }
    }

    #ifdef MPI
    MPI_Allreduce(MPI_IN_PLACE, npart, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, mx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(MPI_IN_PLACE, my, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(MPI_IN_PLACE, mz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
    
    if (*npart>0.0){
        *mx /= *npart;
        *my /= *npart;
        *mz /= *npart;    
        for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
            double *r6 = r_in+c*6;
            if (!atIsNaN(r6[0])) {
                r6[0] -= (*mx-closed_orbit[0])*gx;
                r6[2] -= (*my-closed_orbit[2])*gy;
                r6[5] -= (*mz-closed_orbit[5])*gz;    
            }
        }
    }
    
    
    /*
    printf("%.8f \n", *mz);
    printf("%.4f \n", closed_orbit[5]);
    printf("%d \n", *npart);
    printf("%.2f \n", gz);
    */
    atFree(mx);
    atFree(my);
    atFree(mz);
    atFree(npart);
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double GX, GY, GZ, *closed_orbit;
        GX=atGetDouble(ElemData,"GX"); check_error();
        GY=atGetDouble(ElemData,"GY"); check_error();
        GZ=atGetDouble(ElemData,"GZ"); check_error();
        closed_orbit = atGetDoubleArray(ElemData,"closed_orbit"); check_error();
        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->GX=GX;
        Elem->GY=GY;
        Elem->GZ=GZ;
        Elem->closed_orbit = closed_orbit;
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

        GX=atGetDouble(ElemData,"GX"); check_error();
        GY=atGetDouble(ElemData,"GY"); check_error();
        GZ=atGetDouble(ElemData,"GZ"); check_error();
        closed_orbit = atGetDoubleArray(ElemData,"closed_orbit");check_error();
                
        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->GX=GX;
        Elem->GY=GY;
        Elem->GZ=GZ;
        Elem->closed_orbit = closed_orbit;
        
    if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix: particle array");
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        FeedbackPass(r_in,num_particles,Elem);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(4,1);
        mxSetCell(plhs[0],0,mxCreateString("GX"));
        mxSetCell(plhs[0],1,mxCreateString("GY"));
        mxSetCell(plhs[0],2,mxCreateString("GZ"));
        mxSetCell(plhs[0],3,mxCreateString("closed_orbit"));
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 2 or 0 arguments");
    }
}
#endif

