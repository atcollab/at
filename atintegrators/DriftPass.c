
#include "atelem.c"
#include "atlalib.c"

struct elem 
{
  double Length;
  double *R1;
  double *R2;
  double *T1;
  double *T2;
  double *EApertures;
  double *RApertures;
};

void DriftPass(double *r_in, double le,
	       const double *T1, const double *T2,
	       const double *R1, const double *R2,
	       double *RApertures, double *EApertures,
	       int num_particles)
/* le - physical length
   r_in - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{
  double *r6;
  int c;
  
  #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) default(shared) shared(r_in,num_particles) private(c,r6)
  for (c = 0; c<num_particles; c++) { /*Loop over particles  */
    r6 = r_in+c*6;
    if(!atIsNaN(r6[0])) {
      /*  misalignment at entrance  */
      if (T1) ATaddvv(r6, T1);
      if (R1) ATmultmv(r6, R1);
      /* Check physical apertures at the entrance of the magnet */
      if (RApertures) checkiflostRectangularAp(r6,RApertures);
      if (EApertures) checkiflostEllipticalAp(r6,EApertures);
      ATdrift6(r6, le);
      /* Check physical apertures at the exit of the magnet */
      if (RApertures) checkiflostRectangularAp(r6,RApertures);
      if (EApertures) checkiflostEllipticalAp(r6,EApertures);
      /* Misalignment at exit */
      if (R2) ATmultmv(r6, R2);
      if (T2) ATaddvv(r6, T2);
    }
  }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                double *r_in, int num_particles, struct parameters *Param)
{
/*  if (ElemData) {*/
        if (!Elem) {
            double Length;
            double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
            Length=atGetDouble(ElemData,"Length"); check_error();
            R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
            R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
            T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
            T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
            EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
            RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
            Elem = (struct elem*)atMalloc(sizeof(struct elem));
            Elem->Length=Length;
            Elem->R1=R1;
            Elem->R2=R2;
            Elem->T1=T1;
            Elem->T2=T2;
            Elem->EApertures=EApertures;
            Elem->RApertures=RApertures;
        }
        DriftPass(r_in, Elem->Length, Elem->T1, Elem->T2, Elem->R1, Elem->R2, Elem->RApertures, Elem->EApertures, num_particles);
/*  }
    else {
         atFree(Elem->T1);
         atFree(Elem->T2);
         atFree(Elem->R1);
         atFree(Elem->R2);
         atFree(Elem->EApertures);
         atFree(Elem->RApertures);
     }*/
    return Elem;
}

MODULE_DEF(DriftPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        Length=atGetDouble(ElemData,"Length"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        DriftPass(r_in, Length, T1, T2, R1, R2, RApertures, EApertures, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(1,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(6,1);
            mxSetCell(plhs[1],0,mxCreateString("T1"));
            mxSetCell(plhs[1],1,mxCreateString("T2"));
            mxSetCell(plhs[1],2,mxCreateString("R1"));
            mxSetCell(plhs[1],3,mxCreateString("R2"));
            mxSetCell(plhs[1],4,mxCreateString("RApertures"));
            mxSetCell(plhs[1],5,mxCreateString("EApertures"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
