/* ThinMPolePass
   Accelerator Toolbox 
   A.Terebilo terebilo@ssrl.slac.stanford.edu
*/



#include "mex.h"
#include "elempass.h"
#include "atlalib.c"




void strthinkick(double* r, double* A, double* B,  int max_order)
/******************************************************************************
Calculate and apply a multipole kick to a 6-dimentional
phase space vector in a thin element

IMPORTANT !!! The reference coordinate system can be straight or curved.

If reference coordinate system is straight, the field expansion may still
contain dipole terms: 
PolynomA(1), PolynomB(1) - in MATLAB notation,
A[0], B[0] - C,C++ notation

If reference coordinate system is curved
the dipole field Bo that provides this curvature MUST NOT be included in the dipole term
PolynomB[1](MATLAB notation)(C: B[0] in this function) of the By field expansion
Note: 


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
	double ReSumTemp;
    	for(i=max_order-1;i>=0;i--)
        {   ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
            ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
            ReSum = ReSumTemp;
        }

    r[1] -=  ReSum;
    r[3] += ImSum;
}





void ThinMPolePass(double *r, double *A, double *B, int max_order, 
 double bax, double bay, int num_particles)
{	int c;
	double *r6;   
	for(c = 0;c<num_particles;c++)	/* Loop over particles */
			{		r6 = r+c*6;
			        if(!mxIsNaN(r6[0]))
                    {   strthinkick(r6, A, B, max_order);
                        {   r6[1] += bax*r[4];
                            r6[3] -= bay*r[4];
                            r6[5] -= bax*r[0]-bay*r[2]; /* Path lenghtening */
                        }
                         
                    }                       
			}
}



ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
								double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 4


{	double *A , *B, *tmpdblptr, bax, bay;
	int max_order;
	int *returnptr;
	int *NewFieldNumbers,fnum;
    mxArray *tmpmxptr;

	
	switch(mode)
		{	
		    /* case NO_LOCAL_COPY:	Not used in AT1.3 */
			
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
					
                    /* Optional field 'BendingAngle' */
                    fnum = mxGetFieldNumber(ElemData,"BendingAngle");
	                NewFieldNumbers[3] = fnum;
					if(fnum<0)
                    {   bax = 0;
                        bay = 0;
                    }
					else
                    {   tmpmxptr = mxGetFieldByNumber(ElemData,0,fnum);
                        tmpdblptr = mxGetPr(tmpmxptr);
                        bax = tmpdblptr[0];
                        if(mxGetNumberOfElements(tmpmxptr)>1)
                            bay = tmpdblptr[1];
                        else
                            bay = 0;
                    }
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
				    /* Optional field 'BendingAngle' */
                    if(FieldNumbers[3]<0)
					{   bax = 0;
                        bay = 0;
                    }
					else
					{   tmpmxptr = mxGetFieldByNumber(ElemData,0,3);
                        tmpdblptr = mxGetPr(tmpmxptr);
                        bax = tmpdblptr[0];
                        if(mxGetNumberOfElements(tmpmxptr)>1)
                            bay = tmpdblptr[1];
                        else
                            bay = 0;
                    }
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function ThinMPolePass\n");
				}
		}

	ThinMPolePass(r_in,A, B, max_order,bax,bay,num_particles);
	
	return(returnptr);

}

 






void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	int m,n;
	double *r_in;
	double *A, *B, *tmpdblptr, bax, bay ;  
	int max_order;
	mxArray *tmpmxptr;
    
	if(nrhs)
	{
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
	
    tmpmxptr = mxGetField(prhs[0],0,"BendingAngle");
    
    if(tmpmxptr==NULL)
    {   bax = 0;
        bay = 0;
    }
    else
    {   tmpdblptr = mxGetPr(tmpmxptr);
        bax = tmpdblptr[0];
    if(mxGetNumberOfElements(tmpmxptr)>1)
        bay = tmpdblptr[1];
    else
        bay = 0;
    }

    plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	ThinMPolePass(r_in,A,B,max_order,bax, bay, n);
	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(3,1);
	    mxSetCell(plhs[0],0,mxCreateString("PolynomA"));
	    mxSetCell(plhs[0],1,mxCreateString("PolynomB"));
	    mxSetCell(plhs[0],2,mxCreateString("MaxOrder"));
        if(nlhs>1)  /* optional fields */ 
        {   plhs[1] = mxCreateCellMatrix(1,1);
            mxSetCell(plhs[1],0,mxCreateString("BendingAngle"));
        }
	    
	}

}



