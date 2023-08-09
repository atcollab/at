
#include "atelem.c"
#include "atrandom.c"

struct elem 
{
  double emit_x;
  double emit_y;
  double sigma_dp;
  double tau_x;
  double tau_y;
  double tau_z;
  double beta_x;
  double beta_y;
  double sigma_xp;
  double sigma_yp;
};

void SimpleQuantDiffPass(double *r_in,
           double sigma_xp, double sigma_yp, double sigma_dp,
           double tau_x, double tau_y, double tau_z,
           pcg32_random_t* rng, int num_particles)

{
  double *r6;
  int c, i;
  
  #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) default(shared) shared(r_in,num_particles) private(c,r6)
  for (c = 0; c<num_particles; c++) { /*Loop over particles  */
    r6 = r_in+c*6;
    double randnorm[3];
    for (i = 0; i < 3; i++) {
        randnorm[i] = atrandn_r(rng, 0.0, 1.0);
    }
    
    if(!atIsNaN(r6[0])) {
      if(tau_x>0.0) {
        r6[1] += 2*sigma_xp*sqrt(1/tau_x)*randnorm[0];
      }
      if(tau_y>0.0) {
        r6[3] += 2*sigma_yp*sqrt(1/tau_y)*randnorm[1];
      }
      if(tau_z>0.0) {
        r6[4] += 2*sigma_dp*sqrt(1/tau_z)*randnorm[2];
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
            double emit_x, emit_y, sigma_dp, tau_x, tau_y, tau_z, beta_x, beta_y;
            emit_x=atGetDouble(ElemData,"emit_x"); check_error();
            emit_y=atGetDouble(ElemData,"emit_y"); check_error();
            sigma_dp=atGetDouble(ElemData,"sigma_dp"); check_error();
            tau_x=atGetDouble(ElemData,"tau_x"); check_error();
            tau_y=atGetDouble(ElemData,"tau_y"); check_error();
            tau_z=atGetDouble(ElemData,"tau_z"); check_error();
            beta_x=atGetDouble(ElemData,"beta_x"); check_error();
            beta_y=atGetDouble(ElemData,"beta_y"); check_error();
                            
            Elem = (struct elem*)atMalloc(sizeof(struct elem));
            Elem->emit_x=emit_x;
            Elem->emit_y=emit_y;
            Elem->sigma_dp=sigma_dp;
            Elem->tau_x=tau_x;
            Elem->tau_y=tau_y;
            Elem->tau_z=tau_z;
            Elem->beta_x=beta_x;
            Elem->beta_y=beta_y;
            Elem->sigma_xp=sqrt(emit_x/beta_x);
            Elem->sigma_yp=sqrt(emit_y/beta_y);
            
        }
        SimpleQuantDiffPass(r_in, Elem->sigma_xp, Elem->sigma_yp, Elem->sigma_dp, Elem->tau_x, Elem->tau_y, Elem->tau_z, Param->thread_rng, num_particles);
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

MODULE_DEF(SimpleQuantDiffPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double emit_x, emit_y, sigma_dp, tau_x, tau_y, tau_z, beta_x, beta_y, sigma_xp, sigma_yp;

        emit_x=atGetDouble(ElemData,"emit_x"); check_error();
        emit_y=atGetDouble(ElemData,"emit_y"); check_error();
        sigma_dp=atGetDouble(ElemData,"sigma_dp"); check_error();
        tau_x=atGetDouble(ElemData,"tau_x"); check_error();
        tau_y=atGetDouble(ElemData,"tau_y"); check_error();
        tau_z=atGetDouble(ElemData,"tau_z"); check_error();
        beta_x=atGetDouble(ElemData,"beta_x"); check_error();  
        beta_y=atGetDouble(ElemData,"beta_y"); check_error();  
        sigma_xp=sqrt(emit_x/beta_x);
        sigma_yp=sqrt(emit_y/beta_y);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        SimpleQuantDiffPass(r_in, sigma_xp, sigma_yp, sigma_dp, tau_x, tau_y, tau_z, &pcg32_global, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(8,1);
        mxSetCell(plhs[0],0,mxCreateString("emit_x"));
        mxSetCell(plhs[0],1,mxCreateString("emit_y"));
        mxSetCell(plhs[0],2,mxCreateString("sigma_dp"));
        mxSetCell(plhs[0],3,mxCreateString("tau_x"));
        mxSetCell(plhs[0],4,mxCreateString("tau_y"));
        mxSetCell(plhs[0],5,mxCreateString("tau_z"));
        mxSetCell(plhs[0],6,mxCreateString("beta_x"));
        mxSetCell(plhs[0],7,mxCreateString("beta_y"));

        
        
        /*
        if (nlhs>1) {
            plhs[1] = mxCreateCellMatrix(6,1);
            mxSetCell(plhs[1],0,mxCreateString("T1"));
            mxSetCell(plhs[1],1,mxCreateString("T2"));
            mxSetCell(plhs[1],2,mxCreateString("R1"));
            mxSetCell(plhs[1],3,mxCreateString("R2"));
            mxSetCell(plhs[1],4,mxCreateString("RApertures"));
            mxSetCell(plhs[1],5,mxCreateString("EApertures"));
        }*/
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
