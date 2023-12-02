
#include "atelem.c"
#include "atrandom.c"

struct elem 
{
  double *damp_mat_diag;
  double betax;
  double betay;
  double alphax;
  double alphay;
  double dispx;
  double dispxp;
  double dispy;
  double dispyp;
  double U0;
  double EnergyLossFactor;
};

void SimpleRadiationPass(double *r_in,
           double *damp_mat_diag, double betax, double alphax,
           double betay, double alphay, double dispx, double dispxp,
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
            double U0, EnergyLossFactor, betax, betay, alphax, alphay;
            double dispx, dispxp, dispy, dispyp;
            double *damp_mat_diag;
            
            damp_mat_diag=atGetDoubleArray(ElemData,"damp_mat_diag"); check_error();
            U0=atGetDouble(ElemData,"U0"); check_error();
            alphax=atGetDouble(ElemData,"alphax"); check_error();
            alphay=atGetDouble(ElemData,"alphay"); check_error();
            betax=atGetDouble(ElemData,"betax"); check_error();
            betay=atGetDouble(ElemData,"betay"); check_error();
            dispx=atGetDouble(ElemData,"dispx"); check_error();
            dispxp=atGetDouble(ElemData,"dispxp"); check_error();
            dispy=atGetDouble(ElemData,"dispy"); check_error();
            dispyp=atGetDouble(ElemData,"dispyp"); check_error();
                        
            Elem = (struct elem*)atMalloc(sizeof(struct elem));
            Elem->alphax=alphax;
            Elem->alphay=alphay;
            Elem->betax=betax;
            Elem->betay=betay;
            Elem->U0=U0;
            Elem->dispx=dispx;
            Elem->dispxp=dispxp;
            Elem->dispy=dispy;
            Elem->dispyp=dispyp;
            Elem->EnergyLossFactor=U0/Param->energy;
            Elem->damp_mat_diag=damp_mat_diag;
        }
        SimpleRadiationPass(r_in, Elem->damp_mat_diag, Elem->betax, Elem->alphax, Elem->betay, Elem->alphay, Elem->dispx, Elem->dispxp, Elem->dispy, Elem->dispyp, Elem->EnergyLossFactor, num_particles);
    return Elem;
}

MODULE_DEF(SimpleRadiationPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double alphax, alphay, betax, betay, U0, EnergyLossFactor;
        double dispx, dispxp, dispy, dispyp;
        double *damp_mat_diag;

        damdamp_mat_diagp_mat=atGetDoubleArray(ElemData,"damp_mat_diag"); check_error();
        alphax=atGetDouble(ElemData,"alphax"); check_error();
        alphay=atGetDouble(ElemData,"alphay"); check_error();
        betax=atGetDouble(ElemData,"betax"); check_error();
        betay=atGetDouble(ElemData,"betay"); check_error();
        dispx=atGetDouble(ElemData,"dispx"); check_error();
        dispy=atGetDouble(ElemData,"dispy"); check_error();
        dispxp=atGetDouble(ElemData,"dispxp"); check_error();
        dispyp=atGetDouble(ElemData,"dispyp"); check_error();
        U0=atGetDouble(ElemData,"U0"); check_error();
        EnergyLossFactor=U0/6e9;
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        SimpleRadiationPass(r_in, damp_mat, betax, alphax, betay, alphay, dispx, dispxp, dispy, dispyp, EnergyLossFactor, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(11,1);
        mxSetCell(plhs[0],0,mxCreateString("damp_mat_diag"));
        mxSetCell(plhs[0],1,mxCreateString("alphax"));
        mxSetCell(plhs[0],3,mxCreateString("alphay"));
        mxSetCell(plhs[0],4,mxCreateString("betax"));
        mxSetCell(plhs[0],5,mxCreateString("betay"));
        mxSetCell(plhs[0],6,mxCreateString("dispx"));
        mxSetCell(plhs[0],7,mxCreateString("dispy"));
        mxSetCell(plhs[0],8,mxCreateString("dispxp"));
        mxSetCell(plhs[0],9,mxCreateString("dispyp"));
        mxSetCell(plhs[0],10,mxCreateString("U0"));        
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
