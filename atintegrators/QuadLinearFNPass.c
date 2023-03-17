/* QuadLinearFPass.c 
   Accelerator Toolbox 
   Xiaobiao Huang, SLAC
   Created 3/11/2011
   this is modified from Terebilo's QuadLinearPass to include the quadrupole fringe field effects.
*/


#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "elempass.h"
#include "atlalib.c"


/******************************************************************************/
/* PHYSICS SECTION ************************************************************/

void quad6 (double *r, double L, double K)
{	/* K - is the quadrupole strength defined as
	   (e/Eo)(dBz/dx) [1/m^2] 
	   another notation: g0 [DESY paper]
	*/
	
	double p_norm = 1/(1+r[4]);
	double x, xpr, y ,ypr, g, t ,lt;
	double M12,M21,M34,M43,MVD,MHD;  /* non-0 elements of transfer matrix */
		
	if(K==0) /* Track as a drift */
	{	ATdrift6(r,L);
		return;
	}   
	
	x   = r[0];
	xpr = r[1]*p_norm;
	y   = r[2];
	ypr = r[3]*p_norm;
	
	g  = fabs(K)/(1+r [4]);
	t  = sqrt(g);
	lt = L*t;
				
   if(K>0)
		{	/* Horizontal  */
				MHD = cos(lt); 		
				M12 = sin(lt)/t;
				M21 = -M12*g;		
   			/* Vertical */
				MVD = cosh(lt);		
				M34 = sinh(lt)/t;
				M43 = M34*g;			
		}
	else 
		{	/* Horizontal  */
				MHD = cosh(lt);		
				M12 = sinh(lt)/t;
				M21 = M12*g;			
			/* Vertical */
				MVD = cos(lt); 		
				M34 = sin(lt)/t;
				M43 = -M34*g;			
		}			
	

	/* M transformation compute change in tau first */
	r[0]=  MHD*x + M12*xpr;
	r[1]= (M21*x + MHD*xpr)/p_norm; 
    r[2]=  MVD*y + M34*ypr;
	r[3]= (M43*y + MVD*ypr)/p_norm;  
  
    /* no change in r[4] (delta) */
	
	r[5]+= g*(x*x*(L-MHD*M12)-y*y*(L-MVD*M34))/4;
	r[5]+= (xpr*xpr*(L+MHD*M12)+ypr*ypr*(L+MVD*M34))/4;
	r[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2;
}

/*
Leading term of the linear map for quadrupole fringe field, which is just the phase space scaling,
See Irwin Pac95
*/
void quadlinearfringe(double *r, double I1)
{
	double p_norm = 1/(1+r[4]);
	I1 = I1*p_norm;
	
	r[0] = r[0]*exp(I1); 
	r[1] = r[1]*exp(-I1);
	r[2] = r[2]*exp(-I1);
	r[3] = r[3]*exp(I1);

}

/*
Leading term of the linear map for quadrupole fringe field, which is just the phase space scaling
See Forest & Milutinovic, NIMA 269 (1988) 474
*/

void quadnonlinearfringe(double *r, double kv)
{
	double p_norm = 1/(1+r[4]);
	
	/* R1 reprensets a rotation of the xOy axes by -pi/4 about the s-axis */
	double R1[] ={1/sqrt(2.0), 0, 1/sqrt(2.0), 0, 0, 0,    0, 1/sqrt(2.0),  0, 1/sqrt(2.0), 0, 0,    
		-1/sqrt(2.0), 0, 1/sqrt(2.0), 0, 0, 0,   0, -1/sqrt(2.0), 0, 1/sqrt(2.0),0,0,
	0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0, 1};
	double R2[] ={1/sqrt(2.0), 0, -1/sqrt(2.0), 0, 0, 0,    0, 1/sqrt(2.0),  0, -1/sqrt(2.0), 0, 0,    
		1/sqrt(2.0), 0, 1/sqrt(2.0), 0, 0, 0,   0, 1/sqrt(2.0), 0, 1/sqrt(2.0),0,0,
	  0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0, 1};
	double x, px, y ,py;
			
	
	ATmultmv(r,R2);	 
	x   = r[0];
	px  = r[1];
	y   = r[2];
	py  = r[3];
		 
	r[0] = x - kv/6.*y*y*y*p_norm;
	r[3] = py+ kv/2.*y*y*px*p_norm;
	
	x   = r[0];
	py  = r[3];
	
	r[1] = px + kv/2.*x*x*py*p_norm;
	r[2] = y -  kv/6.*x*x*x*p_norm;
	
  ATmultmv(r,R1);	 
}



void QuadLinearFPass(double *r, double le, double kv, double I1a, double I1b,double *T1, double *T2, double *R1, double *R2, int num_particles)
{	int c;
	double *r6;
	bool useT1, useT2, useR1, useR2;
	char strWarning[256];
	
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
	
	/*sprintf(strWarning, "%f\t%f\t%f\t%f\n",le,kv,I1a,I1b);    
  mexWarnMsgTxt(strWarning);
  */

	for(c = 0;c<num_particles;c++)
		{	r6 = r+c*6;
		    if(!mxIsNaN(r6[0]) & mxIsFinite(r6[4]))
		    /* 
		       function quad6 internally calculates the square root
			   of the energy deviation of the particle 
			   To protect against DOMAIN and OVERFLOW error, check if the
			   fifth component of the phase spacevector r6[4] is finite
			*/
		    {
			/* Misalignment at entrance */
	        if(useT1)
			    ATaddvv(r6,T1);
			if(useR1)
			    ATmultmv(r6,R1);
			
			if (I1a!=0)
				quadnonlinearfringe(r6,kv);
	
			quadlinearfringe(r6, -kv*I1a);
			
			quad6(r6,le,kv);
			
			quadlinearfringe(r6, kv*I1b);
			if (I1b!=0)
				quadnonlinearfringe(r6,-kv);
			
			/* Misalignment at exit */	
			if(useR2)
			    ATmultmv(r6,R2);
		    if(useT2)   
			    ATaddvv(r6,T2);
			}
			
		}		
}

/********** END PHYSICS SECTION ***********************************************/
/******************************************************************************/

/********** WINDOWS DLL GATEWAY SECTION ***************************************/


ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
								double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 8

{	double *pr1, *pr2, *pt1, *pt2 , le, kv, I1a, I1b;   
	int *returnptr,fnum;
	int *NewFieldNumbers;
    switch(mode)
		{	case NO_LOCAL_COPY:	/* Obsolete in AT1.3 et fields by names from MATLAB workspace  */
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
					kv = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[1]));
                    
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
	                
					fnum = mxGetFieldNumber(ElemData,"I1a");
					NewFieldNumbers[6] = fnum;
					if(fnum<0) 
					    I1a = 0;
					else
					   I1a = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[6]));
					
					fnum = mxGetFieldNumber(ElemData,"I1b");
					NewFieldNumbers[7] = fnum;
					if(fnum<0) 
					    I1b = 0;
					else
					   I1b = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[7]));

					returnptr = NewFieldNumbers;

				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
										The second argument ponter to the array of field 
										numbers is previously created with 
										QuadLinPass( ..., MAKE_LOCAL_COPY)
								    */
											
				{	le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					kv = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));

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

					if(FieldNumbers[6]<0)
					    I1a = 0;
					else 
					    I1a = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
					    
					if(FieldNumbers[7]<0)
					    I1b = 0;
					else 
					    I1b = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));    
					                        
					returnptr = FieldNumbers;

		
				}	break;

		}

	QuadLinearFPass(r_in, le, kv, I1a,I1b,pt1, pt2, pr1, pr2 , num_particles);
	return(returnptr);	
}


/********** END WINDOWS DLL GATEWAY SECTION ***************************************/
/********** MATLAB GATEWAY  ***************************************/

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	int m, n;
	double *r_in;   
	double  *pr1, *pr2, *pt1, *pt2 , le, kv, I1a, I1b; 
	mxArray *tmpmxptr;
	

	if(nrhs)
	{
	
	/* ALLOCATE memory for the output array of the same size as the input */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		{mexErrMsgTxt("Second argument must be a 6 x N matrix");}	
	
	/* Required Fields */
	
	tmpmxptr = mxGetField(prhs[0],0,"Length");
	if(tmpmxptr)
	    le = mxGetScalar(tmpmxptr);
	else
	    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
	    
	
	tmpmxptr = mxGetField(prhs[0],0,"K");
	if(tmpmxptr)
	    kv = mxGetScalar(tmpmxptr);
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

  /*I1 for the entrance fringe, 
  I1a and I1b are defined as |\int (K(s)/K0-1)(s-s0) ds|
  value of I1 defined by J. Irwin (PAC95)  would be -K0*I1a for entrance fringe and K0*I1b for exit fringe */
	tmpmxptr = mxGetField(prhs[0],0,"I1a"); 
	if(tmpmxptr)
	    I1a=mxGetScalar(tmpmxptr);
	else
	    I1a=0;
	    
	tmpmxptr = mxGetField(prhs[0],0,"I1b");  /*I1 for the exit fringe,   */
	if(tmpmxptr)
	    I1b=mxGetScalar(tmpmxptr);
	else
	    I1b=0;
	    
	    
	plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	QuadLinearFPass(r_in, le, kv, I1a,I1b,pt1, pt2, pr1, pr2, n);
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