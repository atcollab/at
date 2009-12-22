/* CavityPass.c
   Accelerator Toolbox 
   Revision 3/10/04
   A.Terebilo terebilo@ssrl.slac.stanford.edu
*/

#include "mex.h"
#include "elempass.h"
#include <math.h>
#define TWOPI  6.28318530717959
#define C0  	2.99792458e8 


void CavityPass(double *r_in, double le, double nv, double freq, double lag, int num_particles)
/* le - physical length
   nv - peak voltage (V) normalized to the design enegy (eV)
   r is a 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{	int c, c6;
	double halflength , p_norm, NormL;	
	if(le == 0)
		{	for(c = 0;c<num_particles;c++)
			{	c6 = c*6;
			    if(!mxIsNaN(r_in[c6]))
				    r_in[c6+4] += -nv*sin(TWOPI*freq*(r_in[c6+5]-lag)/C0);
			}
		}
	else
		{	halflength = le/2;
			for(c = 0;c<num_particles;c++)
			{	c6 = c*6;
				if(!mxIsNaN(r_in[c6])) 
				{   p_norm = 1/(1+r_in[c6+4]); 				
				    NormL  = halflength*p_norm;
				    /* Prropagate through a drift equal to half cavity length */
				    r_in[c6+0]+= NormL*r_in[c6+1];
   			        r_in[c6+2]+= NormL*r_in[c6+3];
   			        r_in[c6+5]+= NormL*p_norm*(r_in[c6+1]*r_in[c6+1]+r_in[c6+3]*r_in[c6+3])/2;
				
				    /* Longitudinal momentum kick */
				    r_in[c6+4] += -nv*sin(TWOPI*freq*(r_in[c6+5]-lag)/C0);
				    p_norm = 1/(1+r_in[c6+4]); 				
				    NormL  = halflength*p_norm;
				    /* Prropagate through a drift equal to half cavity length */
				    r_in[c6+0]+= NormL*r_in[c6+1];
   			        r_in[c6+2]+= NormL*r_in[c6+3];
   			        r_in[c6+5]+= NormL*p_norm*(r_in[c6+1]*r_in[c6+1]+r_in[c6+3]*r_in[c6+3])/2;
                }
			}
		}

}




ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
				double *r_in, int num_particles, int mode)


#define NUM_FIELDS_2_REMEMBER 5

{	double le, volt, freq, energy, lag;
	int *returnptr;
	int *NewFieldNumbers, fnum;

    
	switch(mode)
		{	case NO_LOCAL_COPY:	/* Obsolete in AT1.3  */
				{	
				   
				}	break;	
			
			case MAKE_LOCAL_COPY: 	/* Find field numbers first 
									    Save a list of field number in an array
										 and make returnptr point to that array
								    */
				{	
					NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
					
					fnum = mxGetFieldNumber(ElemData,"Length");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;
					
					fnum = mxGetFieldNumber(ElemData,"Voltage");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Voltage' was not found in the element data structure"); 
					NewFieldNumbers[1] = fnum;
					
					fnum = mxGetFieldNumber(ElemData,"Energy");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Energy' was not found in the element data structure"); 
					NewFieldNumbers[2] = fnum;
					
					
					fnum = mxGetFieldNumber(ElemData,"Frequency");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Frequency' was not found in the element data structure"); 
					NewFieldNumbers[3] = fnum;
                    
                    /* Optional field TimeLag */
                    fnum = mxGetFieldNumber(ElemData,"TimeLag");
                    NewFieldNumbers[4] = fnum;
					if(fnum<0) 
					    lag = 0;
                    else
                        lag = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));
					volt = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[1]));
					energy = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[2]));
					freq = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[3]));

					returnptr = NewFieldNumbers;
				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
										 The second argument ponter to the array of field 
										 numbers is previously created with 
										 QuadLinPass( ..., MAKE_LOCAL_COPY)
								    */
											
				{	le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					volt = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
					energy = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
					freq = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
                    /* Optional field TimeLag */
                    if(FieldNumbers[4]<0) 
					    lag = 0;
					else
                        lag = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
                    
					returnptr = FieldNumbers;
				}	break;

			default:
				{	mexErrMsgTxt("No match found for calling mode in function CavityPass\n");
				}
	}

	
	
	CavityPass(r_in, le, volt/energy, freq, lag, num_particles);
	return(returnptr);
}









void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	double volt,freq, energy, lag;
	int m,n;
	double *r_in, le;   
	mxArray *tmpmxptr;
	
	if(nrhs)
	{

	/* ALLOCATE memory for the output array of the same size as the input */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		mexErrMsgTxt("Second argument must be a 6 x N matrix");

	
	tmpmxptr=mxGetField(prhs[0],0,"Length");
	if(tmpmxptr)
		le = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
				
					    
	tmpmxptr=mxGetField(prhs[0],0,"Voltage");
	if(tmpmxptr)
		volt = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Voltage' was not found in the element data structure");
		
    tmpmxptr=mxGetField(prhs[0],0,"Energy");
	if(tmpmxptr)
		energy = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Energy' was not found in the element data structure");
					    
	tmpmxptr=mxGetField(prhs[0],0,"Frequency");
	if(tmpmxptr)
		freq = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Frequency' was not found in the element data structure"); 
    
    tmpmxptr=mxGetField(prhs[0],0,"TimeLag");
	if(tmpmxptr)
		lag = mxGetScalar(tmpmxptr);
	else
		lag = 0; 


    plhs[0] = mxDuplicateArray(prhs[1]);
    r_in = mxGetPr(plhs[0]);
	CavityPass(r_in, le, volt/energy, freq, lag, n);
    }
    else
    {   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(4,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("Voltage"));
	    mxSetCell(plhs[0],2,mxCreateString("Energy"));
	    mxSetCell(plhs[0],3,mxCreateString("Frequency"));
	    if(nlhs>1) /* optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(1,1); 
            mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
	    }
	}

}