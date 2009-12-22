/* CorrectorPass.c 
   Accelerator Toolbox 
   Revision 10/11/01
   A.Terebilo terebilo@ssrl.slac.stanford.edu
   CorrectorPass expects an element to have a fields:
   Length, KickAngle
*/



#include "mex.h"
#include "elempass.h"


void CorrectorPass(double *r_in, double xkick, double ykick, double len,  int num_particles)
/* xkick, ykick - horizontal and vertical kicks in radiand 
   r - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{	int c, c6;
    double NormL, p_norm;
	if(len==0)
	    for(c = 0;c<num_particles;c++)
		{	c6 = c*6;
		    if(!mxIsNaN(r_in[c6]))
		    {   r_in[c6+1] += xkick;
   		        r_in[c6+3] += ykick; 		    
   		    }
		}	
	else
        for(c = 0;c<num_particles;c++)
		{	c6 = c*6;
		    if(!mxIsNaN(r_in[c6]))
		    {   p_norm = 1/(1+r_in[c6+4]); 
			    NormL  = len*p_norm;
	            r_in[c6+5] += NormL*p_norm*(xkick*xkick/3 + ykick*ykick/3 +
   		                r_in[c6+1]*r_in[c6+1] + r_in[c6+3]*r_in[c6+3] + 
   		                r_in[c6+1]*xkick + r_in[c6+3]*ykick)/2;
   		                
			    r_in[c6]   += NormL*(r_in[c6+1]+xkick/2);
		        r_in[c6+1] += xkick;
		        r_in[c6+2] += NormL*(r_in[c6+3]+ykick/2);
   		        r_in[c6+3] += ykick;
   		    }
		}	
}


#define NUM_FIELDS_2_REMEMBER 2

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
				double *r_in, int num_particles, int mode)



{	double *kickptr;
    double le;
	int *returnptr;
	int *NewFieldNumbers, fnum;
    mxArray *tmpmxptr;
    
	switch(mode)
		{	case NO_LOCAL_COPY:	/* Obsolete in AT 1.3  */
				{	returnptr = NULL;
				}	break;	
			
			case MAKE_LOCAL_COPY: 	/* Find field numbers first
										Save a list of field number in an array
										and make returnptr point to that array
									*/
				{	
					NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
					
					fnum = mxGetFieldNumber(ElemData,"KickAngle");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'KickAngle' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;
					
						
					
					fnum = mxGetFieldNumber(ElemData,"Length");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
					NewFieldNumbers[1] = fnum;
										
					
					kickptr = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[1]));
					
					returnptr = NewFieldNumbers;
				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
										The second argument ponter to the array of field 
										numbers is previously created with 
										QuadLinPass( ..., MAKE_LOCAL_COPY)
									*/
											
				{	kickptr = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
				    le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
					returnptr = FieldNumbers;
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function CorrectorPass\n");
				}
	}
	CorrectorPass(r_in, kickptr[0], kickptr[1], le, num_particles);
	return(returnptr);
}


void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	double *kickptr,le;
	int m,n;
	double *r_in;   
	mxArray *tmpmxptr;
	
	if(nrhs)
	{
	/* ALLOCATE memory for the output array of the same size as the input */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6)
		mexErrMsgTxt("Second argument must be a 6 x N matrix");
   
	
	tmpmxptr=mxGetField(prhs[0],0,"KickAngle");
	if(tmpmxptr)
        kickptr = mxGetPr(tmpmxptr);
    else
		mexErrMsgTxt("Required field 'KickAngle' was not found in the element data structure"); 

				
	tmpmxptr=mxGetField(prhs[0],0,"Length");
	if(tmpmxptr)
		le = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 

	
	plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	CorrectorPass(r_in, kickptr[0], kickptr[1],le, n);	
    }
    else
     {   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(2,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("KickAngle"));
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
	    }
	}

}