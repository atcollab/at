/* 
 *  EnergyRampPass.c
 *  Accelerator Toolbox 
 *  25/07/2024
 *  Nicola Carmignani
 */

#include "atconstants.h"
#include "atelem.c"


struct elem
{
    double E0;
    int NPointsRamp;
    double* TurnsRamp;
    double* EnergyRamp;
    /* internal variable */
    int CurrentRampIndex;
};

void EnergyRampPass(double *r_in, struct elem *Elem, int nturn, 
        double *Energy, int num_particles)
{
    double E_entrance=*Energy;
    double E_exit;
    int t1, t2, endramp=0;
    double E1, E2;
    if (Elem->CurrentRampIndex>=Elem->NPointsRamp)
        endramp=1;
    if (nturn >= Elem->TurnsRamp[Elem->CurrentRampIndex] && !endramp)
        Elem->CurrentRampIndex++;
    if (endramp)
        E_exit = Elem->EnergyRamp[Elem->NPointsRamp-1];
    else
    {
        /* linear interpolation */
        if (Elem->CurrentRampIndex==0)
        {
            t1=0;
            E1 = Elem->E0;
        }
        else
        {
            t1 = Elem->TurnsRamp[Elem->CurrentRampIndex-1];
            E1 = Elem->EnergyRamp[Elem->CurrentRampIndex-1];
        }
        t2 = Elem->TurnsRamp[Elem->CurrentRampIndex];
        E2 = Elem->EnergyRamp[Elem->CurrentRampIndex];
        E_exit = E1 + ((E2-E1)/(t2-t1))*(nturn-t1);
    }
    *Energy=E_exit;

    /* loop in the particles to update delta coordinate */
    for (int c = 0; c<num_particles; c++) { /*Loop over particles  */
        /* Get the pointer to the current particle's state vector. */
        double *r6 = r_in + 6*c;
        if (!atIsNaN(r6[0])) {
            double p_norm;
            double xpr, ypr;

            /* calculate angles from tranverse momenta 	*/
            p_norm = 1.0 + r6[4];
            xpr = r6[1]/p_norm;
            ypr = r6[3]/p_norm;
            
            /* change of delta coordinate */
            r6[4] = (E_entrance * (1.0 + r6[4]) - E_exit) / E_exit;
            
            /* recalculate momenta from angles after change of delta */
            p_norm = 1.0 + r6[4];
            r6[1] = xpr*p_norm;
            r6[3] = ypr*p_norm;
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{
    double Energy_entrance=Param->energy;
    double Energy_exit, E0;
    int nturn=Param->nturn;
    
    if (!Elem) {
        int NPointsRamp;
        double *TurnsRamp, *EnergyRamp;
        int CurrentRampIndex=0;
        
        E0 = atGetDouble(ElemData,"E0"); check_error();
        EnergyRamp = atGetDoubleArray(ElemData,"EnergyRamp"); check_error();
        NPointsRamp = atGetLong(ElemData,"NPointsRamp"); check_error();
        TurnsRamp = atGetDoubleArray(ElemData,"TurnsRamp"); check_error();
        
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        
        Elem->NPointsRamp=NPointsRamp;
        Elem->TurnsRamp=TurnsRamp;
        Elem->EnergyRamp=EnergyRamp;
        Elem->CurrentRampIndex=CurrentRampIndex;
    }
    EnergyRampPass(r_in, Elem, nturn, &Param->energy,
            num_particles);
    return Elem;
}

MODULE_DEF(EnergyRampPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	
  if(nrhs == 2)
    {
      double *r_in;
      const mxArray *ElemData = prhs[0];
      int num_particles = mxGetN(prhs[1]);
      double *EnergyRamp, *TurnsRamp, E0;
      int NPointsRamp;
      struct elem *Elem;
      EnergyRamp = atGetDoubleArray(ElemData,"EnergyRamp"); check_error();
      NPointsRamp = atGetLong(ElemData,"NPointsRamp"); check_error();
      E0 = atGetDouble(ElemData,"E0"); check_error();
      TurnsRamp = atGetDoubleArray(ElemData,"TurnsRamp"); check_error();  
      Elem = (struct elem*)atMalloc(sizeof(struct elem));
      
      Elem->NPointsRamp=NPointsRamp;
      Elem->TurnsRamp=TurnsRamp;
      Elem->EnergyRamp=EnergyRamp;
      Elem->CurrentRampIndex=0;
      Elem->E0=E0;
      
      int nturn=0;
      double Energy=200e6;
      
      if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
      /* ALLOCATE memory for the output array of the same size as the input  */
      plhs[0] = mxDuplicateArray(prhs[1]);
      r_in = mxGetDoubles(plhs[0]);
      EnergyRampPass(r_in, Elem, nturn, &Energy, num_particles);
    }
  else if (nrhs == 0)
  {   /* return list of required fields */
      plhs[0] = mxCreateCellMatrix(4,1);
      mxSetCell(plhs[0],0,mxCreateString("E0"));
      mxSetCell(plhs[0],1,mxCreateString("EnergyRamp"));
      mxSetCell(plhs[0],2,mxCreateString("NPointsRamp"));
      mxSetCell(plhs[0],3,mxCreateString("TurnsRamp"));
      if(nlhs>1) /* optional fields */
      {
          plhs[1] = mxCreateCellMatrix(0,1);
      }
  }
  else
  {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
  }
  
}
#endif /* MATLAB_MEX_FILE */
