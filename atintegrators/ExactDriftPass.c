#include "atelem.c"
#include "atlalib.c"
#include "exactdrift.c"

struct elem {
  double Length;
  double *R1;
  double *R2;
  double *T1;
  double *T2;
  double *EApertures;
  double *RApertures;
};

static void drift_pass(double *r_in, double le, const double *T1, const double *T2,
               const double *R1, const double *R2, double *RApertures,
               double *EApertures, int num_particles)
{
  double *r6;
  int c;

  #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD * 10) \
                       default(shared) \
                       shared(r_in, num_particles) \
                       private(c, r6)
  for (c = 0; c < num_particles; c++) { /*Loop over particles  */
    r6 = r_in + c * 6;
    if (!atIsNaN(r6[0])) {

      /*  misalignment at entrance  */
      if (T1) ATaddvv(r6, T1);
      if (R1) ATmultmv(r6, R1);

      /* Check physical apertures at the entrance of the magnet */
      if (RApertures) checkiflostRectangularAp(r6, RApertures);
      if (EApertures) checkiflostEllipticalAp(r6, EApertures);

      exact_drift(r6, le);

      /* Convert absolute path length to path lengthening */
      r6[5] -= le;

      /* Check physical apertures at the exit of the magnet */
      if (RApertures) checkiflostRectangularAp(r6, RApertures);
      if (EApertures) checkiflostEllipticalAp(r6, EApertures);

      /* Misalignment at exit */
      if (R2) ATmultmv(r6, R2);
      if (T2) ATaddvv(r6, T2);
    }
  }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData, struct elem *Elem,
                                      double *r_in, int num_particles,
                                      struct parameters *Param) {
  if (!Elem) {
    double Length = atGetDouble(ElemData, "Length"); check_error();
    double *R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
    double *R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
    double *T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
    double *T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();
    double *EApertures = atGetOptionalDoubleArray(ElemData, "EApertures"); check_error();
    double *RApertures = atGetOptionalDoubleArray(ElemData, "RApertures"); check_error();
    Elem = (struct elem *)atMalloc(sizeof(struct elem));
    Elem->Length = Length;
    Elem->R1 = R1;
    Elem->R2 = R2;
    Elem->T1 = T1;
    Elem->T2 = T2;
    Elem->EApertures = EApertures;
    Elem->RApertures = RApertures;
  }
  drift_pass(r_in, Elem->Length, Elem->T1, Elem->T2, Elem->R1, Elem->R2,
             Elem->RApertures, Elem->EApertures, num_particles);
  return Elem;
}

MODULE_DEF(ExactDriftPass) /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs >= 2) {
    double *r_in;
    const mxArray *ElemData = prhs[0];
    int num_particles = mxGetN(prhs[1]);

    double Length = atGetDouble(ElemData, "Length"); check_error();
    double *R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
    double *R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
    double *T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
    double *T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();
    double *EApertures = atGetOptionalDoubleArray(ElemData, "EApertures"); check_error();
    double *RApertures = atGetOptionalDoubleArray(ElemData, "RApertures"); check_error();
    if (mxGetM(prhs[1]) != 6)
      mexErrMsgIdAndTxt("AT:WrongArg",
                        "Second argument must be a 6 x N matrix");
    /* ALLOCATE memory for the output array of the same size as the input  */
    plhs[0] = mxDuplicateArray(prhs[1]);
    r_in = mxGetDoubles(plhs[0]);
    drift_pass(r_in, Length, T1, T2, R1, R2, RApertures, EApertures,
               num_particles);
  } else if (nrhs == 0) {
    /* list of required fields */
    int i0 = 0;
    plhs[0] = mxCreateCellMatrix(1, 1);
    mxSetCell(plhs[0], i0++, mxCreateString("Length"));
    if (nlhs > 1) {
      /* list of optional fields */
      int i1 = 0;
      plhs[1] = mxCreateCellMatrix(6, 1);
      mxSetCell(plhs[1], i1++, mxCreateString("T1"));
      mxSetCell(plhs[1], i1++, mxCreateString("T2"));
      mxSetCell(plhs[1], i1++, mxCreateString("R1"));
      mxSetCell(plhs[1], i1++, mxCreateString("R2"));
      mxSetCell(plhs[1], i1++, mxCreateString("RApertures"));
      mxSetCell(plhs[1], i1++, mxCreateString("EApertures"));
    }
  } else {
    mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
  }
}
#endif /*defined(MATLAB_MEX_FILE)*/
