/* 
 *  RFCavityPass.c
 *  Accelerator Toolbox 
 *  22/09/2015
 *  Nicola Carmignani
 */

#include "atconstants.h"
#include "attrackfunc.c"


struct elem 
{
  double Length;
  double Voltage;
  double Energy;
  double Frequency;
  double HarmNumber;
  double TimeLag;
  double PhaseLag;
  /* variables for voltage ramp */
  int NPointsRamp;
  double* TurnsRamp;
  double* VoltageRamp;
  /* internal variable */
  int CurrentRampIndex;
};

void RFCavityPass(double *r_in, double le, double Voltage, double energy, 
        double freq, double h, double lag, double philag, int nturn, 
        double T0, int NPointsRamp, double *TurnsRamp, double *VoltageRamp, 
        int* CurrentRampIndex, int num_particles)
/* le - physical length
   nv - peak voltage (V) normalized to the design enegy (eV)
   r is a 6-by-N matrix of initial conditions reshaped into
   1-d array of 6*N elements
*/
{
    double nv;
    if (NPointsRamp)
    {
        int  endramp=0;
        double t1, t2, V1, V2;
        if (*CurrentRampIndex>=NPointsRamp)
            endramp=1;
        if (nturn >= TurnsRamp[*CurrentRampIndex] && !endramp)
            *CurrentRampIndex=*CurrentRampIndex+1;
        if (endramp)
            Voltage = VoltageRamp[NPointsRamp-1];
        else
        {
            /* linear interpolation */
            if (*CurrentRampIndex==0)
            {
                t1=0;
                V1 = Voltage;
            }
            else
            {
                t1 = TurnsRamp[*CurrentRampIndex-1];
                V1 = VoltageRamp[*CurrentRampIndex-1];
            }
            t2 = TurnsRamp[*CurrentRampIndex];
            V2 = VoltageRamp[*CurrentRampIndex];
            Voltage = V1 + ((V2-V1)/(t2-t1))*(nturn-t1);
        }
    }
    
    nv=Voltage/energy;
    trackRFCavity(r_in, le, nv, freq, h, lag, philag, nturn, T0, num_particles);
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{
    int nturn=Param->nturn;
    double T0=Param->T0;
    if (!Elem) {
        double Length, Voltage, Energy, Frequency, TimeLag, PhaseLag, *VoltageRamp, *TurnsRamp;
        int NPointsRamp;
        Length=atGetDouble(ElemData,"Length"); check_error();
        Voltage=atGetDouble(ElemData,"Voltage"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        PhaseLag=atGetOptionalDouble(ElemData,"PhaseLag",0); check_error();
        VoltageRamp = atGetOptionalDoubleArray(ElemData,"VoltageRamp"); check_error();
        NPointsRamp = atGetOptionalLong(ElemData,"NPointsRamp",0); check_error();
        TurnsRamp = atGetOptionalDoubleArray(ElemData,"TurnsRamp"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->Voltage=Voltage;
        Elem->Energy=Param->energy;
        Elem->Frequency=Frequency;
        Elem->HarmNumber=round(Frequency*T0);
        Elem->TimeLag=TimeLag;
        Elem->PhaseLag=PhaseLag;
        Elem->NPointsRamp=NPointsRamp;
        Elem->TurnsRamp=TurnsRamp;
        Elem->VoltageRamp=VoltageRamp;
        Elem->CurrentRampIndex=0;
    }
    
    RFCavityPass(r_in, Elem->Length, Elem->Voltage, Param->energy, Elem->Frequency, 
            Elem->HarmNumber, Elem->TimeLag, Elem->PhaseLag, nturn, T0, 
            Elem->NPointsRamp, Elem->TurnsRamp, Elem->VoltageRamp, 
            &Elem->CurrentRampIndex,num_particles);
    return Elem;
}

MODULE_DEF(RFCavityPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	
  if(nrhs == 2)
    {
      double *r_in;
      const mxArray *ElemData = prhs[0];
      int num_particles = mxGetN(prhs[1]);
      double Length=atGetDouble(ElemData,"Length");
      double Voltage=atGetDouble(ElemData,"Voltage");
      double Energy=atGetDouble(ElemData,"Energy");
      double Frequency=atGetDouble(ElemData,"Frequency");
      double TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0);
      double PhaseLag=atGetOptionalDouble(ElemData,"PhaseLag",0);
      double T0=1.0/Frequency;      /* Does not matter since nturns == 0 */
      double HarmNumber=round(Frequency*T0);
      double *VoltageRamp = atGetOptionalDoubleArray(ElemData,"VoltageRamp");
      int NPointsRamp = atGetOptionalLong(ElemData,"NPointsRamp",0);
      double *TurnsRamp = atGetOptionalDoubleArray(ElemData,"TurnsRamp");
      if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
      /* ALLOCATE memory for the output array of the same size as the input  */
      plhs[0] = mxDuplicateArray(prhs[1]);
      r_in = mxGetDoubles(plhs[0]);
      RFCavityPass(r_in, Length, Voltage, Energy, Frequency, HarmNumber,
              TimeLag, PhaseLag, 0, T0, NPointsRamp, TurnsRamp, VoltageRamp,
              0, num_particles);
    }
  else if (nrhs == 0)
  {   /* return list of required fields */
      plhs[0] = mxCreateCellMatrix(4,1);
      mxSetCell(plhs[0],0,mxCreateString("Length"));
      mxSetCell(plhs[0],1,mxCreateString("Voltage"));
      mxSetCell(plhs[0],2,mxCreateString("Energy"));
      mxSetCell(plhs[0],3,mxCreateString("Frequency"));
      if(nlhs>1) /* optional fields */
      {
          plhs[1] = mxCreateCellMatrix(5,1);
          mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
          mxSetCell(plhs[1],1,mxCreateString("PhaseLag"));
          mxSetCell(plhs[1],2,mxCreateString("NPointsRamp"));
          mxSetCell(plhs[1],3,mxCreateString("TurnsRamp"));
          mxSetCell(plhs[1],4,mxCreateString("VoltageRamp"));
      }
  }
  else
  {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
  }
  
}
#endif /* MATLAB_MEX_FILE */
