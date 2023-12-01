
#include "atelem.c"
#include "atrandom.c"

struct elem 
{
  double *damp_mat;
  double betax;
  double betay;
  double alphax;
  double alphay;
  double dispx;
  double dispy;
  double U0;
  double EnergyLossFactor;
};

void SimpleRadiationPass(double *r_in,
           double *damp_mat, double betax, double alphax,
           double betay, double alphay, double dispx,
           double dispy, double EnergyLossFactor, int num_particles)

{
  double *r6;
  int c, i;
  double *dmatx = damp_mat;
  double *dmatxp = damp_mat + 6;
  double *dmaty = damp_mat + 2*6;
  double *dmatyp = damp_mat + 3*6;
  double *dmatdp = damp_mat + 4*6;
  double *dmatz = damp_mat + 5*6;
  double x,xp,y,yp,z,dpp;
      
  #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) default(shared) shared(r_in,num_particles) private(c,r6)
  for (c = 0; c<num_particles; c++) { /*Loop over particles  */
    r6 = r_in+c*6;
    
    if(!atIsNaN(r6[0])) {

      dpp = r6[4];
      x = r6[0]/sqrt(betax);
      xp = r6[1]/(1.0+dpp);
      y = r6[2]/sqrt(betay);
      yp = r6[3]/(1.0+dpp);
      z = r6[5];        
    
      if(dmatx[0]!=1.0) {
        x *= dmatx[0];
      }

      if(dmatxp[1]!=1.0) {
        xp *= dmatxp[1];
      }
      
      if(dmaty[2]!=1.0) {
        y *= dmaty[2];
      }

      if(dmatyp[3]!=1.0) {
        yp *= dmatyp[3];
      }
      
      if(dmatdp[4]!=1.0) {
        dpp *= dmatdp[4];
      }

      if(dmatz[5]!=1.0) {
        z *= dmatz[5];
      }

      
      r6[0] = x*sqrt(betax);
      r6[1] = xp*(1+dpp);
      r6[2] = y*sqrt(betay);
      r6[3] = yp*(1+dpp);
      r6[4] = dpp - EnergyLossFactor;
      r6[5] = z;
      
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
            double dispx, dispy;
            double *damp_mat;
            
            damp_mat=atGetDoubleArray(ElemData,"damp_mat"); check_error();
            U0=atGetDouble(ElemData,"U0"); check_error();
            alphax=atGetDouble(ElemData,"alphax"); check_error();
            alphay=atGetDouble(ElemData,"alphay"); check_error();
            betax=atGetDouble(ElemData,"betax"); check_error();
            betay=atGetDouble(ElemData,"betay"); check_error();
            dispx=atGetDouble(ElemData,"dispx"); check_error();
            dispy=atGetDouble(ElemData,"dispy"); check_error();
            
            Elem = (struct elem*)atMalloc(sizeof(struct elem));
            Elem->alphax=alphax;
            Elem->alphay=alphay;
            Elem->betax=betax;
            Elem->betay=betay;
            Elem->U0=U0;
            Elem->dispx=dispx;
            Elem->dispy=dispy;
            Elem->EnergyLossFactor=U0/Param->energy;
            Elem->damp_mat=damp_mat;
        }
        SimpleRadiationPass(r_in, Elem->damp_mat, Elem->betax, Elem->alphax, Elem->betay, Elem->alphay, Elem->dispx, Elem->dispy, Elem->EnergyLossFactor, num_particles);
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
        double *damp_mat;

        damp_mat=atGetDoubleArray(ElemData,"damp_mat"); check_error();
        alphax=atGetDouble(ElemData,"alphax"); check_error();
        alphay=atGetDouble(ElemData,"alphay"); check_error();
        betax=atGetDouble(ElemData,"betax"); check_error();
        betay=atGetDouble(ElemData,"betay"); check_error();
        dispx=atGetDouble(ElemData,"dispx"); check_error();
        dispy=atGetDouble(ElemData,"dispy"); check_error();
        U0=atGetDouble(ElemData,"U0"); check_error();
        EnergyLossFactor=U0/6e9;
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        SimpleRadiationPass(r_in, damp_mat, betax, alphax, betay, alphay, dispx, dispy, EnergyLossFactor, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(9,1);
        mxSetCell(plhs[0],0,mxCreateString("damp_mat"));
        mxSetCell(plhs[0],1,mxCreateString("alphax"));
        mxSetCell(plhs[0],3,mxCreateString("alphay"));
        mxSetCell(plhs[0],4,mxCreateString("betax"));
        mxSetCell(plhs[0],5,mxCreateString("betay"));
        mxSetCell(plhs[0],6,mxCreateString("dispx"));
        mxSetCell(plhs[0],7,mxCreateString("dispy"));
        mxSetCell(plhs[0],8,mxCreateString("U0"));        
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
