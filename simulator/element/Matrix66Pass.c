/* DriftPass.c 
   Accelerator Toolbox 
   Revision 6/26/00
   A.Terebilo terebilo@ssrl.slac.stanford.edu
*/

#include "mex.h"
#include "elempass.h"
#include "atlalib.c"

void Matrix66Pass(double *r, double *M, int num_particles)

{	int c, c6;
    
	for(c = 0;c<num_particles;c++)
	{   c6 = c*6;
	    if(!mxIsNaN(r[c6]))
	        ATmultmv(r+c*6,M);
	}

}




ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
				double *r_in, int num_particles, int mode)


#define NUM_FIELDS_2_REMEMBER 1

{	double *M;

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
					fnum = mxGetFieldNumber(ElemData,"M66");
					if(fnum<0)
					    mexErrMsgTxt("Required field 'M66' was not found in the element data structure"); 
					else
					{   NewFieldNumbers[0] = fnum;
					    M = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					    returnptr = NewFieldNumbers;
					}
				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
									   The second argument ponter to the array of field 
									    numbers is previously created with 
									     QuadLinPass( ..., MAKE_LOCAL_COPY)
									*/
											
				{	M = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					returnptr = FieldNumbers;
				}	break;
	}
	Matrix66Pass(r_in, M, num_particles);
	return(returnptr);
}









void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	int m,n;
	double *r_in, *M;
	mxArray *tmpmxptr;
		
	if(nrhs)
	{
	/* ALLOCATE memory for the output array of the same size as the input  */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		mexErrMsgTxt("Second argument must be a 6 x N matrix");
    
    tmpmxptr=mxGetField(prhs[0],0,"M66");
	if(tmpmxptr)
	    M = mxGetPr(tmpmxptr);
    else
	    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
		
		
    plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	
	Matrix66Pass(r_in,M, n);	
	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(1,1);
	    mxSetCell(plhs[0],0,mxCreateString("M66"));
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
	    }
	}

}
