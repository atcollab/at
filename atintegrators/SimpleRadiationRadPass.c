
#include "atelem.c"
#include "atrandom.c"

struct elem 
{
  double *damp_mat_diag;
  double dispx;
  double dispxp;
  double dispy;
  double dispyp;
  double U0;
  double EnergyLossFactor;
};

void SimpleRadiationRadPass(double *r_in,
           double *damp_mat_diag, double dispx, double dispxp,
           double dispy, double dispyp, double EnergyLossFactor, int num_particles)

{
  double *r6;
  int c, i;
  double x,xp,y,yp,z,dpp;
  
  #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) default(shared) shared(r_in,num_particles) private(c,r6)
  for (c = 0; c<num_particles; c++) { /*Loop over particles  */
    r6 = r_in+c*6;
    
    if(!atIsNaN(r6[0])) {

      dpp = r6[4];
      x = r6[0] - dispx*dpp;
      xp = r6[1] - dispxp*dpp;
      
      y = r6[2] - dispy*dpp;
      yp = r6[3] - dispyp*dpp;
      
      z = r6[5];        
    
      if(damp_mat_diag[0]!=1.0) {
        x *= damp_mat_diag[0];
      }

      if(damp_mat_diag[1]!=1.0) {
        xp *= damp_mat_diag[1];
      }
      
      if(damp_mat_diag[2]!=1.0) {
        y *= damp_mat_diag[2];
      }

      if(damp_mat_diag[3]!=1.0) {
        yp *= damp_mat_diag[3];
      }
      
      if(damp_mat_diag[4]!=1.0) {
        dpp *= damp_mat_diag[4];
      }

      if(damp_mat_diag[5]!=1.0) {
        z *= damp_mat_diag[5];
      }

      r6[4] = dpp;
      r6[5] = z;
      r6[0] = x + dispx*dpp;
      r6[1] = xp + dispxp*dpp;
      r6[2] = y + dispy*dpp;
      r6[3] = yp + dispyp*dpp;      
      
      r6[4] -= EnergyLossFactor;      
    }
  }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                double *r_in, int num_particles, struct parameters *Param)
{
/*  if (ElemData) {*/
        if (!Elem) {
            double U0, EnergyLossFactor;
            double dispx, dispxp, dispy, dispyp;
            double *damp_mat_diag;
            
            damp_mat_diag=atGetDoubleArray(ElemData,"damp_mat_diag"); check_error();
            U0=atGetDouble(ElemData,"U0"); check_error();
            dispx=atGetOptionalDouble(ElemData,"dispx",0.0); check_error();
            dispxp=atGetOptionalDouble(ElemData,"dispxp",0.0); check_error();
            dispy=atGetOptionalDouble(ElemData,"dispy",0.0); check_error();
            dispyp=atGetOptionalDouble(ElemData,"dispyp",0.0); check_error();
                        
            Elem = (struct elem*)atMalloc(sizeof(struct elem));
            Elem->U0=U0;
            Elem->dispx=dispx;
            Elem->dispxp=dispxp;
            Elem->dispy=dispy;
            Elem->dispyp=dispyp;
            Elem->EnergyLossFactor=U0/Param->energy;
            Elem->damp_mat_diag=damp_mat_diag;
        }
        SimpleRadiationRadPass(r_in, Elem->damp_mat_diag, Elem->dispx, Elem->dispxp, Elem->dispy, Elem->dispyp, Elem->EnergyLossFactor, num_particles);
    return Elem;
}

MODULE_DEF(SimpleRadiationRadPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double U0, EnergyLossFactor;
        double dispx, dispxp, dispy, dispyp;
        double *damp_mat_diag;

        damp_mat_diag=atGetDoubleArray(ElemData,"damp_mat_diag"); check_error();
        dispx=atGetOptionalDouble(ElemData,"dispx",0.0); check_error();
        dispy=atGetOptionalDouble(ElemData,"dispy",0.0); check_error();
        dispxp=atGetOptionalDouble(ElemData,"dispxp",0.0); check_error();
        dispyp=atGetOptionalDouble(ElemData,"dispyp",0.0); check_error();
        U0=atGetDouble(ElemData,"U0"); check_error();
        EnergyLossFactor=U0/6e9;
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        SimpleRadiationRadPass(r_in, damp_mat_diag, dispx, dispxp, dispy, dispyp, EnergyLossFactor, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(2,1);
        mxSetCell(plhs[0],0,mxCreateString("damp_mat_diag"));
        mxSetCell(plhs[0],1,mxCreateString("U0"));
        
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(4,1); /* No optional fields */
            mxSetCell(plhs[1],0,mxCreateString("dispx"));
            mxSetCell(plhs[1],1,mxCreateString("dispxp"));
            mxSetCell(plhs[1],2,mxCreateString("dispy"));
            mxSetCell(plhs[1],3,mxCreateString("dispyp"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
