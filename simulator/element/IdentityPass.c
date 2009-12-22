/* IdentityPass.c 
   Accelerator Toolbox
   Revision 7/16/03
   A.Terebilo terebilo@slac.stanford.edu
*/

#include "mex.h"
#include "elempass.h"


/* void IdentityPass(double *r_in)
{	
}
*/




ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
				double *r_in, int num_particles, int mode)


{	/* IdentityPass(r_in); */
	return(NULL);
}


void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	int m,n;
	if(nrhs)
	{
	/* ALLOCATE memory for the output array of the same size as the input */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		mexErrMsgTxt("Second argument must be a 6 x N matrix");
   
	plhs[0] = mxDuplicateArray(prhs[1]);

	/*
	r_in = mxGetPr(plhs[0]);
	IdentityPass(r_in);	
    */
    }
    else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(0,0);
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
	    }
	}
}