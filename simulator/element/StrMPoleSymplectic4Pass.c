#include "mex.h"
#include "elempass.h"
#include "atlalib.c"


#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656



void strthinkick(double* r, double* A, double* B, double L, int max_order)
/***************************************************************************** 
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
	double ReSumTemp;
    	for(i=max_order-1;i>=0;i--)
        {   ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
            ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
            ReSum = ReSumTemp;
        }

    r[1] -=  L*ReSum;
    r[3] +=  L*ImSum;
}

void fastdrift(double* r, double NormL)

/*   NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
     in the loop if momentum deviation (delta) does not change
     such as in 4-th order symplectic integrator w/o radiation
*/

{   double dx = NormL*r[1];
    double dy = NormL*r[3];
    r[0]+= dx;
    r[2]+= dy;
    r[5]+= NormL*(r[1]*r[1]+r[3]*r[3])/(2*(1+r[4]));
}

void StrMPoleSymplectic4Pass(double *r, double le, double *A, double *B,
					int max_order, int num_int_steps,
					double *T1, double *T2,	
					double *R1, double *R2, int num_particles)
{	int c,m;
	double norm, NormL1, NormL2;	
	double *r6;
	bool useT1, useT2, useR1, useR2;
	double SL, L1, L2, K1, K2;
	SL = le/num_int_steps;
	L1 = SL*DRIFT1;
	L2 = SL*DRIFT2;
	K1 = SL*KICK1;
	K2 = SL*KICK2;
	
	if(T1==NULL)
	    useT1=false;
	else 
	    useT1=true;  
	    
    if(T2==NULL)
	    useT2=false; 
	else 
	    useT2=true;  
	
	if(R1==NULL)
	    useR1=false; 
	else 
	    useR1=true;  
	    
    if(R2==NULL)
	    useR2=false;
	else 
	    useR2=true;
	for(c = 0;c<num_particles;c++)	/*Loop over particles  */
			{	r6 = r+c*6;	
			    if(!mxIsNaN(r6[0]))
			    {   
					/*  misalignment at entrance  */
					if(useT1)
			            ATaddvv(r6,T1);
			        if(useR1)
			            ATmultmv(r6,R1);
					/*  integrator  */
					for(m=0; m < num_int_steps; m++) /*  Loop over slices */
						{	r6 = r+c*6;	
							norm = 1/(1+r6[4]);
							NormL1 = L1*norm;
							NormL2 = L2*norm;
							fastdrift(r6, NormL1);
           					strthinkick(r6, A, B,  K1, max_order);
           					fastdrift(r6, NormL2);
           					strthinkick(r6, A, B, K2, max_order);
           					fastdrift(r6, NormL2);
           					strthinkick(r6, A, B,  K1, max_order);
           					fastdrift(r6, NormL1);	
						}  
					
					/* Misalignment at exit */	
			        if(useR2)
			            ATmultmv(r6,R2);
		            if(useT2)   
			            ATaddvv(r6,T2);
			    }
			}
}



ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
								double *r_in, int num_particles, int mode)
#define NUM_FIELDS_2_REMEMBER 9
{	int fnum;
	double *A , *B;
	double  *pr1, *pr2, *pt1, *pt2;   
	int max_order, num_int_steps;
	double le;
	int *returnptr;
	int *NewFieldNumbers;

	switch(mode)
		{	case NO_LOCAL_COPY:	/* NOT used in AT1.3 Get fields by names from MATLAB workspace   */
				{	
				
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

					/*  Populate */
					
					fnum = mxGetFieldNumber(ElemData,"PolynomA");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'PolynomA' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;
					A = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					fnum = mxGetFieldNumber(ElemData,"PolynomB");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'PolynomB' was not found in the element data structure"); 
					NewFieldNumbers[1] = fnum;
					B = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					
					fnum = mxGetFieldNumber(ElemData,"MaxOrder");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'MaxOrder' was not found in the element data structure"); 
					NewFieldNumbers[2] = fnum;
					max_order = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					fnum = mxGetFieldNumber(ElemData,"NumIntSteps");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'NumIntSteps' was not found in the element data structure"); 
					NewFieldNumbers[3] = fnum;
					num_int_steps = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					fnum = mxGetFieldNumber(ElemData,"Length");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
					NewFieldNumbers[4] = fnum;
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					fnum = mxGetFieldNumber(ElemData,"R1");
					NewFieldNumbers[5] = fnum;
					if(fnum<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					

					fnum = mxGetFieldNumber(ElemData,"R2");
					NewFieldNumbers[6] = fnum;
					if(fnum<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
                    fnum = mxGetFieldNumber(ElemData,"T1");
	                NewFieldNumbers[7] = fnum;
					if(fnum<0)
					    pt1 = NULL;
					else
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
	                
	                fnum = mxGetFieldNumber(ElemData,"T2");
	                NewFieldNumbers[8] = fnum;
					if(fnum<0)
					    pt2 = NULL;
					else
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
								
					
								
				
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
					num_int_steps = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
					
					/* Optional fields */
					if(FieldNumbers[5]<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
					
					if(FieldNumbers[6]<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
					
					    
					if(FieldNumbers[7]<0)
					    pt1 = NULL;
					else    
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));
					    
					if(FieldNumbers[8]<0)
					    pt2 = NULL;
					else 
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[8]));
					returnptr = FieldNumbers;
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function StrMPoleSymplectic4Pass\n");
				}
		}

    StrMPoleSymplectic4Pass(r_in, le, A, B, max_order, num_int_steps,
									pt1, pt2, pr1, pr2, num_particles);
	return(returnptr);
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	int m,n;
	double *r_in;
	double le, *A, *B, *pr1, *pr2, *pt1, *pt2;  
	int max_order, num_int_steps;
	mxArray *tmpmxptr;

	if(nrhs)
	{
    /* ALLOCATE memory for the output array of the same size as the input  */
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
				        
	tmpmxptr = mxGetField(prhs[0],0,"NumIntSteps");
	if(tmpmxptr)   
		num_int_steps = (int)mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'NumIntSteps' was not found in the element data structure");    
				    
	tmpmxptr = mxGetField(prhs[0],0,"Length");
	if(tmpmxptr)
		le = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Length' was not found in the element data structure");    
					
				    
	/* Optionnal arguments */    
	tmpmxptr = mxGetField(prhs[0],0,"R1");
	if(tmpmxptr)
	    pr1 = mxGetPr(tmpmxptr);
	else
	    pr1=NULL; 
	    
	tmpmxptr = mxGetField(prhs[0],0,"R2");
	if(tmpmxptr)
	    pr2 = mxGetPr(tmpmxptr);
	else
	    pr2=NULL; 

	tmpmxptr = mxGetField(prhs[0],0,"T1");

	
	if(tmpmxptr)
	    pt1=mxGetPr(tmpmxptr);
	else
	    pt1=NULL;

	tmpmxptr = mxGetField(prhs[0],0,"T2");
	if(tmpmxptr)
	    pt2=mxGetPr(tmpmxptr);
	else
	    pt2=NULL;
	
	
	plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	StrMPoleSymplectic4Pass(r_in, le, A, B, max_order, num_int_steps,
									pt1, pt2, pr1, pr2, n);
	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(5,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("PolynomA"));
	    mxSetCell(plhs[0],2,mxCreateString("PolynomB"));
	    mxSetCell(plhs[0],3,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],4,mxCreateString("NumIntSteps"));	    
	    
        if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(4,1);
	        mxSetCell(plhs[1],0,mxCreateString("T1"));
	        mxSetCell(plhs[1],1,mxCreateString("T2"));
	        mxSetCell(plhs[1],2,mxCreateString("R1"));
	        mxSetCell(plhs[1],3,mxCreateString("R2"));
	    }
    }
}



