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
};

void RFCavityPass(double *r_in, double le, double nv, double freq, double h, double lag, double philag,
                  int nturn, double T0, int num_particles)
/* le - physical length
   nv - peak voltage (V) normalized to the design enegy (eV)
   r is a 6-by-N matrix of initial conditions reshaped into
   1-d array of 6*N elements
*/
{
    trackRFCavity(r_in, le, nv, freq, h, lag, philag, nturn, T0, num_particles);
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{
    int nturn=Param->nturn;
    double T0=Param->T0;
    double energy = Param->energy;
    if (!Elem) {
        double Length, Voltage, Energy, Frequency, TimeLag, PhaseLag;
        Length=atGetDouble(ElemData,"Length"); check_error();
        Voltage=atGetDouble(ElemData,"Voltage"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        /*optional fields*/
        Energy=atGetOptionalDouble(ElemData,"Energy",energy); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        PhaseLag=atGetOptionalDouble(ElemData,"PhaseLag",0); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->Voltage=Voltage;
        Elem->Energy=Energy;
        Elem->Frequency=Frequency;
        Elem->HarmNumber=round(Frequency*T0);
        Elem->TimeLag=TimeLag;
        Elem->PhaseLag=PhaseLag;
    }
    if (energy == 0.0) energy = Elem->Energy;

    RFCavityPass(r_in, Elem->Length, Elem->Voltage/energy, Elem->Frequency, Elem->HarmNumber, Elem->TimeLag,
                 Elem->PhaseLag, nturn, T0, num_particles);
    return Elem;
}

MODULE_DEF(RFCavityPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	
  if (nrhs >= 2) {
      double Energy = 0.0;
      double rest_energy = 0.0;
      double charge = -1.0;
      double *r_in;
      const mxArray *ElemData = prhs[0];
      int num_particles = mxGetN(prhs[1]);
      double Length=atGetDouble(ElemData,"Length");
      double Voltage=atGetDouble(ElemData,"Voltage");
      Energy=atGetOptionalDouble(ElemData,"Energy",Energy);
      double Frequency=atGetDouble(ElemData,"Frequency");
      double TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0);
      double PhaseLag=atGetOptionalDouble(ElemData,"PhaseLag",0);
      double T0=1.0/Frequency;      /* Does not matter since nturns == 0 */
      double HarmNumber=round(Frequency*T0);
      if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);

      if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
      /* ALLOCATE memory for the output array of the same size as the input  */
      plhs[0] = mxDuplicateArray(prhs[1]);
      r_in = mxGetDoubles(plhs[0]);

      RFCavityPass(r_in, Length, Voltage/Energy, Frequency, HarmNumber, TimeLag, PhaseLag, 0, T0, num_particles);

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
          plhs[1] = mxCreateCellMatrix(2,1);
          mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
          mxSetCell(plhs[1],1,mxCreateString("PhaseLag"));
      }
  }
  else
  {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
  }
  
}
#endif /* MATLAB_MEX_FILE */
