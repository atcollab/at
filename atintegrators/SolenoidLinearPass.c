/* SolenoidLinearPass.c 
   Accelerator Toolbox 
   Revision 7/17/03
   A.Terebilo terebilo@slac.stanford.edu
*/

#include "mex.h"
#include "atlalib.c"
#include "elempass.h"
#include <math.h>


void SolenoidLinearPass(double *r_in, double le, double ks, double *T1, double *T2, double *R1, double *R2, int num_particles)
/* Constant field hard edge model is assumed
   le - physical length
   ks - solenoid strength Bo / (B*rho)
   r_in - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{	int c;
	double *r6, p_norm, H, S, C, x, xpr, y, ypr, NormL;
	
	bool useT1, useT2, useR1, useR2;
	
	if (T1==NULL)
	    useT1=false;
	else 
	    useT1=true;  
	    
    if (T2==NULL)
	    useT2=false; 
	else 
	    useT2=true;  
	
	if (R1==NULL)
	    useR1=false;
	else 
	    useR1=true;  
	    
    if (R2==NULL)
	    useR2=false;
	else 
	    useR2=true;
	
	if(ks!=0)
	    for(c = 0;c<num_particles;c++)
		{	r6= r_in+c*6;
			p_norm = 1/(1+r6[4]); 
			
			
			/* Misalignment at entrance */
	        if (useT1)
			    ATaddvv(r6,T1);
			if (useR1)
			    ATmultmv(r6,R1);
			    
			x   = r6[0];
	        xpr = r6[1]*p_norm;
	        y   = r6[2];
	        ypr = r6[3]*p_norm;
			H   = ks*p_norm/2;
			S  = sin(le*H);
			C  = cos(le*H);
			
			
			r6[0]=     x*C*C         +  xpr*C*S/H       +   y*C*S         +   ypr*S*S/H;
	        r6[1]= ( - x*H*C*S       +  xpr*C*C         -   y*H*S*S       +   ypr*C*S         )/p_norm; 
            r6[2]=   - x*C*S         -  xpr*S*S/H       +   y*C*C         +   ypr*C*S/H;
	        r6[3]= (   x*H*S*S       -  xpr*C*S         -   y*C*S*H       +   ypr*C*C         )/p_norm;
   			r6[5]+= le*(H*H*(x*x+y*y) + 2*H*(xpr*y-ypr*x) +xpr*xpr+ypr*ypr)/2;
			
   			/* Misalignment at exit */	
			if (useR2)
			    ATmultmv(r6,R2);
		    if (useT2)   
			    ATaddvv(r6,T2);
   			
   			
		}
    else /* Drift */
        for(c = 0;c<num_particles;c++)
		{	r6= r_in+c*6;
			p_norm = 1/(1+r6[4]);  
			NormL  = le*p_norm;
   			r6[0]+= NormL*r6[1];
   			r6[2]+= NormL*r6[3];
   			r6[5]+= NormL*p_norm*(r6[1]*r6[1]+r6[3]*r6[3])/2;
			
		}
}



ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
				double *r_in, int num_particles, int mode)


#define NUM_FIELDS_2_REMEMBER 6

{	double le, ks, *pr1, *pr2, *pt1, *pt2 ;
	int *returnptr;
	int *NewFieldNumbers, fnum;
    
	
	
	switch(mode)
		{	case NO_LOCAL_COPY:	/* Not used since AT1.3 Get fields by names from MATLAB workspace  */
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
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));
					
					fnum = mxGetFieldNumber(ElemData,"K");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'K' was not found in the element data structure"); 
					NewFieldNumbers[1] = fnum;
					ks = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[1]));
					
					fnum = mxGetFieldNumber(ElemData,"R1");
					NewFieldNumbers[2] = fnum;
					if(fnum<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					

					fnum = mxGetFieldNumber(ElemData,"R2");
					NewFieldNumbers[3] = fnum;
					if(fnum<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
                    fnum = mxGetFieldNumber(ElemData,"T1");
	                NewFieldNumbers[4] = fnum;
					if(fnum<0)
					    pt1 = NULL;
					else
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
	                
	                fnum = mxGetFieldNumber(ElemData,"T2");
	                NewFieldNumbers[5] = fnum;
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
											
				{	le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					ks = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
					/* Optional fields */
					if(FieldNumbers[2]<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
					
					if(FieldNumbers[3]<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
					
					    
					if(FieldNumbers[4]<0)
					    pt1 = NULL;
					else    
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
					    
					if(FieldNumbers[5]<0)
					    pt2 = NULL;
					else 
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
					returnptr = FieldNumbers;
				}	break;

	}

	
	
	SolenoidLinearPass(r_in, le, ks, pt1, pt2, pr1, pr2, num_particles);
	return(returnptr);
}









void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	int m,n;
	double *r_in, le, ks, *pr1, *pr2, *pt1, *pt2;   
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
				
					    
	tmpmxptr=mxGetField(prhs[0],0,"K");
	if(tmpmxptr)
		ks = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'K' was not found in the element data structure"); 
					    
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
	SolenoidLinearPass(r_in, le, ks, pt1, pt2, pr1, pr2, n);
	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(2,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("K"));
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(4,1);
	        mxSetCell(plhs[1],0,mxCreateString("T1"));
	        mxSetCell(plhs[1],1,mxCreateString("T2"));
	        mxSetCell(plhs[1],2,mxCreateString("R1"));
	        mxSetCell(plhs[1],3,mxCreateString("R2"));
	    }
	}


}
