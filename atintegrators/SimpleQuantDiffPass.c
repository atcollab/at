
#include "atelem.c"
#include "atrandom.c"

struct elem 
{
  double emitx;
  double emity;
  double espread;
  double taux;
  double tauy;
  double tauz;
  double betax;
  double betay;
  double sigma_xp;
  double sigma_yp;
  double U0;
  double EnergyLossFactor;
};

void SimpleQuantDiffPass(double *r_in,
           double sigma_xp, double sigma_yp, double espread,
           double taux, double tauy, double tauz, double EnergyLossFactor,
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
      if(taux!=0.0) {
        r6[1] -= 2*r6[1]/taux;
      }
      if(sigma_xp!=0.0) {
        r6[1] += 2*sigma_xp*sqrt(1/taux)*randnorm[0];
      }
      if(sigma_yp!=0.0) {
        r6[3] += 2*sigma_yp*sqrt(1/tauy)*randnorm[1];
      }
      if(tauy!=0.0) {
        r6[3] -= 2*r6[3]/tauy;
      }
      if(espread!=0.0) {
        r6[4] += 2*espread*sqrt(1/tauz)*randnorm[2];
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
            double emitx, emity, espread, taux, tauy, tauz, betax, betay, U0, EnergyLossFactor;
            emitx=atGetDouble(ElemData,"emitx"); check_error();
            emity=atGetDouble(ElemData,"emity"); check_error();
            espread=atGetDouble(ElemData,"espread"); check_error();
            taux=atGetDouble(ElemData,"taux"); check_error();
            tauy=atGetDouble(ElemData,"tauy"); check_error();
            tauz=atGetDouble(ElemData,"tauz"); check_error();
            betax=atGetDouble(ElemData,"betax"); check_error();
            betay=atGetDouble(ElemData,"betay"); check_error();
            U0=atGetDouble(ElemData,"U0"); check_error();
            
            Elem = (struct elem*)atMalloc(sizeof(struct elem));
            Elem->emitx=emitx;
            Elem->emity=emity;
            Elem->espread=espread;
            Elem->taux=taux;
            Elem->tauy=tauy;
            Elem->tauz=tauz;
            Elem->betax=betax;
            Elem->betay=betay;
            Elem->sigma_xp=sqrt(emitx/betax);
            Elem->sigma_yp=sqrt(emity/betay);
            Elem->U0=U0;
            Elem->EnergyLossFactor=U0/Param->energy;
        }
        SimpleQuantDiffPass(r_in, Elem->sigma_xp, Elem->sigma_yp, Elem->espread, Elem->taux, Elem->tauy, Elem->tauz, Elem->EnergyLossFactor, Param->thread_rng, num_particles);
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
        double emitx, emity, espread, taux, tauy, tauz, betax, betay, sigma_xp, sigma_yp, U0, EnergyLossFactor;

        emitx=atGetDouble(ElemData,"emitx"); check_error();
        emity=atGetDouble(ElemData,"emity"); check_error();
        espread=atGetDouble(ElemData,"espread"); check_error();
        taux=atGetDouble(ElemData,"taux"); check_error();
        tauy=atGetDouble(ElemData,"tauy"); check_error();
        tauz=atGetDouble(ElemData,"tauz"); check_error();
        betax=atGetDouble(ElemData,"betax"); check_error();  
        betay=atGetDouble(ElemData,"betay"); check_error();  
        U0=atGetDouble(ElemData,"U0"); check_error();  
        sigma_xp=sqrt(emitx/betax);
        sigma_yp=sqrt(emity/betay);
        EnergyLossFactor=U0/6e9;
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        SimpleQuantDiffPass(r_in, sigma_xp, sigma_yp, espread, taux, tauy, tauz, EnergyLossFactor, &pcg32_global, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(9,1);
        mxSetCell(plhs[0],0,mxCreateString("emitx"));
        mxSetCell(plhs[0],1,mxCreateString("emity"));
        mxSetCell(plhs[0],2,mxCreateString("espread"));
        mxSetCell(plhs[0],3,mxCreateString("taux"));
        mxSetCell(plhs[0],4,mxCreateString("tauy"));
        mxSetCell(plhs[0],5,mxCreateString("tauz"));
        mxSetCell(plhs[0],6,mxCreateString("betax"));
        mxSetCell(plhs[0],7,mxCreateString("betay"));
        mxSetCell(plhs[0],8,mxCreateString("U0"));
        
        
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
