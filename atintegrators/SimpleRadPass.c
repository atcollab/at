
#include "atelem.c"
#include "atrandom.c"

struct elem 
{
  double taux;
  double tauy;
  double tauz;
  double U0;
  double EnergyLossFactor;
};

void SimpleRadPass(double *r_in,
           double taux, double tauy, double tauz,
           double EnergyLossFactor, int num_particles)

{
  double *r6;
  int c, i;
  
  #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) default(shared) shared(r_in,num_particles) private(c,r6)
  for (c = 0; c<num_particles; c++) { /*Loop over particles  */
    r6 = r_in+c*6;
    
    if(!atIsNaN(r6[0])) {
      if(taux!=0.0) {
        r6[1] -= 2*r6[1]/taux;
      }
      if(tauy!=0.0) {
        r6[3] -= 2*r6[3]/tauy;
      }
      if(tauz!=0.0) {
        r6[4] -= 2*r6[4]/tauz;
      }
      if(EnergyLossFactor>0.0) {
        r6[4] -= EnergyLossFactor;
      }
    }
  }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                double *r_in, int num_particles, struct parameters *Param)
{
/*  if (ElemData) {*/
        if (!Elem) {
            double taux, tauy, tauz, U0, EnergyLossFactor;
            taux=atGetDouble(ElemData,"taux"); check_error();
            tauy=atGetDouble(ElemData,"tauy"); check_error();
            tauz=atGetDouble(ElemData,"tauz"); check_error();
            U0=atGetDouble(ElemData,"U0"); check_error();
            
            Elem = (struct elem*)atMalloc(sizeof(struct elem));
            Elem->taux=taux;
            Elem->tauy=tauy;
            Elem->tauz=tauz;
            Elem->U0=U0;
            Elem->EnergyLossFactor=U0/Param->energy;
        }
        SimpleQuantDiffPass(r_in, Elem->taux, Elem->tauy, Elem->tauz, Elem->EnergyLossFactor, num_particles);
    return Elem;
}

MODULE_DEF(SimpleQuantDiffPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double taux, tauy, tauz, U0, EnergyLossFactor;

        taux=atGetDouble(ElemData,"taux"); check_error();
        tauy=atGetDouble(ElemData,"tauy"); check_error();
        tauz=atGetDouble(ElemData,"tauz"); check_error();
        U0=atGetDouble(ElemData,"U0"); check_error();
        EnergyLossFactor=U0/6e9;
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        SimpleQuantDiffPass(r_in, taux, tauy, tauz, EnergyLossFactor, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(4,1);
        mxSetCell(plhs[0],0,mxCreateString("taux"));
        mxSetCell(plhs[0],1,mxCreateString("tauy"));
        mxSetCell(plhs[0],2,mxCreateString("tauz"));
        mxSetCell(plhs[0],3,mxCreateString("U0"));        
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
