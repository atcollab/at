/* DriftPass.c 
   Accelerator Toolbox 
   Revision 6/26/00
   A.Terebilo terebilo@ssrl.slac.stanford.edu
*/

#include "mex.h"
#include "elempass.h"


void DriftPass(double *r_in, double le, int num_particles)
/* le - physical length
   r_in - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{	int c, c6;
	double p_norm, NormL;
	for(c = 0;c<num_particles;c++)
		{	c6 = c*6;
		    if(!mxIsNaN(r_in[c6]))
			{    p_norm = 1/(1+r_in[c6+4]); 
			    NormL  = le*p_norm;
   			    r_in[c6+0]+= NormL*r_in[c6+1];
   			    r_in[c6+2]+= NormL*r_in[c6+3];
   			    r_in[c6+5]+= NormL*p_norm*(r_in[c6+1]*r_in[c6+1]+r_in[c6+3]*r_in[c6+3])/2;
   			 }
			
		}
}




ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
				double *r_in, int num_particles, int mode)


#define NUM_FIELDS_2_REMEMBER 1

{	double le;

	int *returnptr;
	int fnum,*NewFieldNumbers;

	switch(mode)
		{	case NO_LOCAL_COPY:	/* Not used in AT1.3 Get fields by names from MATLAB workspace  */
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
					else
					{   NewFieldNumbers[0] = fnum;
					    le = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));
					    returnptr = NewFieldNumbers;
					}
				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
									   The second argument ponter to the array of field 
									    numbers is previously created with 
									     QuadLinPass( ..., MAKE_LOCAL_COPY)
									*/
											
				{	le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					returnptr = FieldNumbers;
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function DriftPass\n");
				}
	}
	DriftPass(r_in, le, num_particles);
	return(returnptr);
}









void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	double le;
	int m,n;
	double *r_in;
	mxArray *tmpmxptr;
		
	if(nrhs)
	{
	/* ALLOCATE memory for the output array of the same size as the input  */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		mexErrMsgTxt("Second argument must be a 6 x N matrix");
    
    tmpmxptr=mxGetField(prhs[0],0,"Length");
	if(tmpmxptr)
	    le = mxGetScalar(tmpmxptr);
    else
	    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
		
		
    plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	
	DriftPass(r_in,le, n);	
	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(1,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
	    }
	}

}
