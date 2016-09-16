/* 
 *  RFCavityPass.c
 *  Accelerator Toolbox 
 *  22/09/2015
 *  Nicola Carmignani
 */

#include "at.h"
#define TWOPI  6.28318530717959
#define C0  	2.99792458e8 

void RFCavityPass(double *r_in, double le, double nv, double freq, double h, double lag, int nturn, double T0, int num_particles)
/* le - physical length
   nv - peak voltage (V) normalized to the design enegy (eV)
   r is a 6-by-N matrix of initial conditions reshaped into
   1-d array of 6*N elements
*/
{	int c, c6;
  double halflength , p_norm, NormL;
  /* I get T0 and nturn from matlab global variables */
  /*T0=getT0FromMatlab();*/
  
  if(le == 0)
    {
      for(c = 0;c<num_particles;c++)
	{	c6 = c*6;
	  if(!mxIsNaN(r_in[c6]))
	    r_in[c6+4] += -nv*sin(TWOPI*freq*((r_in[c6+5]-lag)/C0 - (h/freq-T0)*nturn ));
	}
    }
  else
    {	halflength = le/2;
      for(c = 0;c<num_particles;c++)
	{	c6 = c*6;
	  if(!mxIsNaN(r_in[c6])) 
	    {   p_norm = 1/(1+r_in[c6+4]); 				
	      NormL  = halflength*p_norm;
	      /* Propagate through a drift equal to half cavity length */
	      r_in[c6+0]+= NormL*r_in[c6+1];
	      r_in[c6+2]+= NormL*r_in[c6+3];
	      r_in[c6+5]+= NormL*p_norm*(r_in[c6+1]*r_in[c6+1]+r_in[c6+3]*r_in[c6+3])/2;
	      /* Longitudinal momentum kick */
	      r_in[c6+4] += -nv*sin(TWOPI*freq*((r_in[c6+5]-lag)/C0 - (h/freq-T0)*nturn ));
	      p_norm = 1/(1+r_in[c6+4]); 				
	      NormL  = halflength*p_norm;
	      /* Propagate through a drift equal to half cavity length */
	      r_in[c6+0]+= NormL*r_in[c6+1];
	      r_in[c6+2]+= NormL*r_in[c6+3];
	      r_in[c6+5]+= NormL*p_norm*(r_in[c6+1]*r_in[c6+1]+r_in[c6+3]*r_in[c6+3])/2;
	    }
	}
    }
} 

#ifdef MATLAB_MEX_FILE

#include "elempass.h"
#include "mxutils.c"

ExportMode int* trackFunction(const mxArray *ElemData,int *FieldNumbers,
			      double *r_in, int num_particles,struct parameters *Param)
#define NUM_FIELDS_2_REMEMBER 6
{	double le, volt, freq, T0, h, energy, lag;
  int nturn;
  nturn=Param->nturn;
  T0=Param->T0;
  /*printf("turn=%d\nT0=%f\nRingLength=%f\n~~~~~\n",Param->nturn,Param->T0,Param->RingLength);*/
  switch(Param->mode)
    {	
    case NO_LOCAL_COPY:	/* Obsolete in AT1.3 */
      {	
      }	break;	
    case MAKE_LOCAL_COPY:	/* Find field numbers first
				   Save a list of field number in an array
				   and make returnptr point to that array
				*/
      /* Allocate memory for integer array of
	 field numbers for faster futurereference
      */
      FieldNumbers = (int *) mxCalloc(NUM_FIELDS_2_REMEMBER, sizeof(int));
      /*  Populate */
      FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "Length");
      FieldNumbers[1] = GetRequiredFieldNumber(ElemData, "Voltage");
      FieldNumbers[2] = GetRequiredFieldNumber(ElemData, "Energy");
      FieldNumbers[3] = GetRequiredFieldNumber(ElemData, "Frequency");
      FieldNumbers[4] = GetRequiredFieldNumber(ElemData, "HarmNumber");
      FieldNumbers[5] = mxGetFieldNumber(ElemData,"TimeLag");
      
      /* Fall through next section... */
    case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
				   The second argument pointer to the array of field
				   numbers is previously created with
				   RFCavityPass(..., MAKE_LOCAL_COPY)
				*/
      le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
      volt = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
      energy = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
      freq = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
      h = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
      /* Optional field TimeLag */
      if(FieldNumbers[5]<0) 
	lag = 0;
      else
	lag = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
      break;
      
    default:
      mexErrMsgTxt("No match found for calling mode in function RFCavityPass\n");
    }
  RFCavityPass(r_in, le, volt/energy, freq, h, lag, nturn, T0, num_particles);
  return FieldNumbers;
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	
  if(nrhs == 2)
    {
      double *r_in;
      mxArray *tmpmxptr;
      double T0;
      double le = mxGetScalar(GetRequiredField(prhs[0], "Length"));	
      int num_particles = mxGetN(prhs[1]);
      double volt = mxGetScalar(GetRequiredField(prhs[0], "Voltage"));	
      double energy = mxGetScalar(GetRequiredField(prhs[0], "Energy"));	
      double freq = mxGetScalar(GetRequiredField(prhs[0], "Frequency"));	
      double h = mxGetScalar(GetRequiredField(prhs[0], "HarmNumber"));	
      /* Optional arguments */
      double lag;
      tmpmxptr=mxGetField(prhs[0],0,"TimeLag");
      if(tmpmxptr)
	lag = mxGetScalar(tmpmxptr);
      else
	lag = 0;
      plhs[0] = mxDuplicateArray(prhs[1]);
      r_in = mxGetPr(plhs[0]);
      T0=h/freq;
      RFCavityPass(r_in, le, volt/energy, freq, h, lag,0, T0, num_particles);
    }
  else if (nrhs == 0) 
    {   /* return list of required fields */
      plhs[0] = mxCreateCellMatrix(5,1);
      mxSetCell(plhs[0],0,mxCreateString("Length"));
      mxSetCell(plhs[0],1,mxCreateString("Voltage"));
      mxSetCell(plhs[0],2,mxCreateString("Energy"));
      mxSetCell(plhs[0],3,mxCreateString("Frequency"));
      mxSetCell(plhs[0],4,mxCreateString("HarmNumber"));
      if(nlhs>1) /* optional fields */ 
	{   
	  plhs[1] = mxCreateCellMatrix(1,1); 
	  mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
	}
    }
  else 
    {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
  
}
#endif /* MATLAB_MEX_FILE */
