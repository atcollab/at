#include "mex.h"
#include <math.h>
#include "elempass.h"


#define SQR(X) ((X)*(X))

void set2zero(double *r6)
{	int i;
	
    for(i=1;i<6;i++)
		r6[i] = 0;
}

void markaslost(double *r6)
{	int i;
    
    r6[0] = mxGetNaN();
	
    for(i=1;i<6;i++)
		r6[i] =0;
}

void AperturePass(double *r_in, double *limitsptr, int num_particles)
{   /*  Checks X and Y of each input 6-vector and marks the corresponding element in 
    lossflag array with 0 if X,Y are exceed the limits given by limitsptr array
	limitsptr has 4 elements: (MinX, MaxX, MinY, MaxY) */
	 
	int i, c, c6;
	for(c = 0;c<num_particles;c++)
	{   c6 = c*6;
        if(!mxIsNaN(r_in[c6])) /*  check if this particle is already marked as lost			*/
        {   /* check limits for X position */
			if(r_in[c6+0]<limitsptr[0] || r_in[c6+0]>limitsptr[1] || r_in[c6+2]<limitsptr[2] || r_in[c6+2]>limitsptr[3])
			    markaslost(r_in+c6);
			else
			    for(i=0;i<6;i++)
                {    if(!mxIsFinite(r_in[c6+i]))
                        {   markaslost(r_in+c6);
                            break;
                        }
                }
         }
    }
}

#ifndef NOMEX

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
				double *r_in, int num_particles, int mode)


#define NUM_FIELDS_2_REMEMBER 1

{	int *returnptr;
	int *NewFieldNumbers, fnum;

	double * limitsptr;

	switch(mode)
		{	case NO_LOCAL_COPY:	/* OBSOLETE ON at 1.3  */
		
				{	
				}	break;	
			
			case MAKE_LOCAL_COPY: 	/* Find field numbers first
									   Save a list of field number in an array
									   and make returnptr point to that array
									*/
				{	
					NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
					fnum = mxGetFieldNumber(ElemData,"Limits");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Limits' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;
					limitsptr = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));

					returnptr = NewFieldNumbers;
				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
										The second argument ponter to the array of field 
										numbers is previously created with 
										QuadLinPass( ..., MAKE_LOCAL_COPY)
									*/	
				{	limitsptr = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					returnptr = FieldNumbers;
				}	break;

			default:
				{	mexErrMsgTxt("No match for calling mode in function AperturePass\n");
				}
	}


	AperturePass(r_in,limitsptr,num_particles);

	return(returnptr);
}













void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	
    double *r_in, *LIMITS;   
  	int m,n;
	mxArray *tmpmxptr;

	if(nrhs)
	{
    /* ALLOCATE memory for the output array of the same size as the input  */
	r_in = mxGetPr(prhs[1]);
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		{mexErrMsgTxt("Second argument must be a 6 x N matrix");}	
	
		


	tmpmxptr = mxGetField(prhs[0],0,"Limits");
	if(tmpmxptr)
	    LIMITS = mxGetPr(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Limits' was not found in the element data structure"); 

    plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);		
		
	AperturePass(r_in,LIMITS,n);

	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(1,1);
	    mxSetCell(plhs[0],0,mxCreateString("Limits"));
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
	    }
	}

}
#endif /*NOMEX*/