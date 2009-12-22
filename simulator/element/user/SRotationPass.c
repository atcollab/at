/* SolenoidPass.c 
   Accelerator Toolbox 
   Revision 7/16/03
   Christoph Steier, CSteier@lbl.gov
*/

#include "mex.h"
#include "elempass.h"
#include "atlalib.c"
#include <stdlib.h>
#include <math.h>

/******************************************************************************/
/* PHYSICS SECTION ************************************************************/

void SRotationPass(double *r, double *R, int num_particles)
{	int c;
	double *r6;
	for(c = 0;c<num_particles;c++)
		{	r6 = r+c*6;	
			
			ATmultmv(r6,R);
						
		}		
}

/********** END PHYSICS SECTION ***********************************************/
/******************************************************************************/

/********** WINDOWS DLL GATEWAY SECTION ***************************************/


ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
								double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 1

{	double *pr;   
    mxArray *tmpmxptr;
	int *returnptr,fnum;
	int *NewFieldNumbers;

	switch(mode)
		{	case NO_LOCAL_COPY:	/* Get fields by names from MATLAB workspace  */
				{						
					tmpmxptr=mxGetField(ElemData,0,"R");
				    if(tmpmxptr)
                        pr = mxGetPr(tmpmxptr);
					else
					    mexErrMsgTxt("Required field 'R' was not found in the element data structure"); 
				
					returnptr = NULL;
				
				}	break;	
			case MAKE_LOCAL_COPY: 	/* Find field numbers first
										Save a list of field number in an array
										and make returnptr point to that array
								    */
				{						
					NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
					
					fnum = mxGetFieldNumber(ElemData,"R");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'R' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;	                

			
					pr = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));

					returnptr = NewFieldNumbers;

				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
										The second argument ponter to the array of field 
										numbers is previously created with 
										SRotationPass( ..., MAKE_LOCAL_COPY)
								    */
											
				{	pr = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));

					returnptr = FieldNumbers;

		
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function SRotationPass\n");
				}
		}

	SRotationPass(r_in,  pr, num_particles);
	return(returnptr);	
}


/********** END WINDOWS DLL GATEWAY SECTION ***************************************/
/********** MATLAB GATEWAY  ***************************************/

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	int m, n;
	double *r_in;   
	double  *pr; 
	mxArray *tmpmxptr;

	/* ALLOCATE memory for the output array of the same size as the input */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		{mexErrMsgTxt("Second argument must be a 6 x N matrix");}	
		
		
	tmpmxptr = mxGetField(prhs[0],0,"R");
	if(tmpmxptr)
	    pr = mxGetPr(tmpmxptr);
	else
	    mexErrMsgTxt("Required field 'R' was not found in the element data structure"); 
	
	plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	SRotationPass(r_in,  pr, n);
}


