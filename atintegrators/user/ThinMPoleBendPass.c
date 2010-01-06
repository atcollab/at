/* ThinMPolePassBend
   Accelerator Toolbox 
   
   Revision 3/10/04
   A.Terebilo terebilo@ssrl.slac.stanford.edu
      
   Created 07/29/2003
   C. Steier CSteier@lbl.gov
   
   interprets dipole term in PolynomB as bending angle instead of
   corrector kick (reference system changes direction).
*/



#include "mex.h"
#include "elempass.h"
#include "atlalib.c"


#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656



void strthinkick(double* r, double* A, double* B,  int max_order)
/******************************************************************************
Calculate and apply a multipole kick to a 6-dimentional
phase space vector in a straight element ( quadrupole)

IMPORTANT !!!
The reference coordinate system is straight but the field expansion may still
contain dipole terms: PolynomA(1), PolynomB(1) - in MATLAB notation,
A[0], B[0] - C,C++ notation


   Note: in the US convention the transverse multipole field is written as:

                         max_order+1
                           ----
                           \                       n-1
	   (B + iB  )/ B rho  =  >   (ia  + b ) (x + iy)
         y    x            /       n    n
	                        ----
                          n=1
	is a polynomial in (x,y) with the highest order = MaxOrder
	

	Using different index notation 
   
                         max_order
                           ----
                           \                       n
	   (B + iB  )/ B rho  =  >   (iA  + B ) (x + iy)
         y    x            /       n    n
	                       ----
                          n=0

	A,B: i=0 ... max_order
   [0] - dipole, [1] - quadrupole, [2] - sextupole ...
   units for A,B[i] = 1/[m]^(i+1)
	Coeficients are stroed in the PolynomA, PolynomB field of the element
	structure in MATLAB

	A[i] (C++,C) =  PolynomA(i+1) (MATLAB) 
	B[i] (C++,C) =  PolynomB(i+1) (MATLAB) 
	i = 0 .. MaxOrder

******************************************************************************/
{  int i;
	double ReSum = B[max_order];
 	double ImSum = A[max_order];
 	double B0temp,A0temp;
	double ReSumTemp;
	
	B0temp=B[0];
	B[0]=0.0;
	A0temp=A[0];
	A[0]=0.0;
	
	/* thin multipole kick: quadrupole and higher */
	
    	for(i=max_order-1;i>=0;i--)
        {   ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
            ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
            ReSum = ReSumTemp;
        }

    r[1] -=  ReSum;
    r[3] += ImSum;
    
    B[0]=B0temp;
    A[0]=A0temp;
    
    /* dispersion prime from energy dependend dipole kick */
    r[1] += B0temp*r[4];
    r[3] -= A0temp*r[4];
    
    /* path lengthening */
    r[5] -= B0temp*r[0]-A0temp*r[2];
}





void ThinMPolePass(double *r, double *A, double *B, int max_order, int num_particles)


{	int c;
	double *r6;   
	
	for(c = 0;c<num_particles;c++)	/* Loop over particles */
			{		r6 = r+c*6;	
			            strthinkick(r6, A, B, max_order);
			}
}



ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
								double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 3


{	double *A , *B;
	int max_order;
	int *returnptr;
	int *NewFieldNumbers,fnum;
	mxArray *tmpmxptr;

	
	switch(mode)
		{	case NO_LOCAL_COPY:	/* Get fields by names from MATLAB workspace  */
				{	tmpmxptr =mxGetField(ElemData,0,"PolynomA");
				    if(tmpmxptr)
				        A = mxGetPr(tmpmxptr);
				    else
				        mexErrMsgTxt("Required field 'PolynomA' was not found in the element data structure"); 
				    
					tmpmxptr =mxGetField(ElemData,0,"PolynomB");
				    if(tmpmxptr)   
					    B = mxGetPr(tmpmxptr);
					else
					    mexErrMsgTxt("Required field 'PolynomB' was not found in the element data structure");
					    
				    tmpmxptr = mxGetField(ElemData,0,"MaxOrder");
				    if(tmpmxptr)
				        max_order = (int)mxGetScalar(tmpmxptr);
					else
					    mexErrMsgTxt("Required field 'MaxOrder' was not found in the element data structure");
					returnptr = NULL;
				
				}	break;	
	
			case MAKE_LOCAL_COPY: 	/* Find field numbers first
									   Save a list of field number in an array
									   and make returnptr point to that array
									*/
				{	
					/* Allocate memory for integer array of 
					  field numbers for faster futurereference
					*/
		
					NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));

					/* Populate */
					fnum = mxGetFieldNumber(ElemData,"PolynomA");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'PolynomA' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;
					
					fnum = mxGetFieldNumber(ElemData,"PolynomB");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'PolynomB' was not found in the element data structure"); 
					NewFieldNumbers[1] = fnum;
					
					
					
					fnum = mxGetFieldNumber(ElemData,"MaxOrder");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'MaxOrder' was not found in the element data structure"); 
					NewFieldNumbers[2] = fnum;
										
					A = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));
					B = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[1]));
					max_order = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[2]));
									
					returnptr = NewFieldNumbers;

				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
									   The second argument ponter to the array of field 
									   numbers is previously created with 
									   QuadLinPass( ..., MAKE_LOCAL_COPY)
									*/
											
				{	A = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					B = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
					max_order = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
				
					returnptr = FieldNumbers;
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function ThinMPolePass\n");
				}
		}

	ThinMPolePass(r_in,A, B, max_order,num_particles);
	
	return(returnptr);

}

 






void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	int m,n;
	double *r_in;
	double *A, *B ;  
	int max_order;
	mxArray *tmpmxptr;


/* ALLOCATE memory for the output array of the same size as the input */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		mexErrMsgTxt("Second argument must be a 6 x N matrix");
	

	tmpmxptr =mxGetField(prhs[0],0,"PolynomA");
	if(tmpmxptr)
	    A = mxGetPr(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'PolynomA' was not found in the element data structure"); 
				    
	tmpmxptr =mxGetField(prhs[0],0,"PolynomB");
	if(tmpmxptr)   
	    B = mxGetPr(tmpmxptr);
	else
	    mexErrMsgTxt("Required field 'PolynomB' was not found in the element data structure");
					    
	tmpmxptr = mxGetField(prhs[0],0,"MaxOrder");
	if(tmpmxptr)
	    max_order = (int)mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'MaxOrder' was not found in the element data structure");
		
    plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);

	ThinMPolePass(r_in,A,B,max_order, n);
	

}



