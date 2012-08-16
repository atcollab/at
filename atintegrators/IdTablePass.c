/* IdTablePass.c 
   Accelerator Toolbox 
   Created: 13/11/08
   Z.Martí zeus@cells.es

   Based in the matlab routine:
	WigTablePass.m - The tracking table is described in
	P. Elleaume, "A new approach to the electron beam dynamics in undulators
	and wigglers", EPAC92.
 
*/

#include <math.h>
#include "mex.h"
#include "elempass.h"
#include "atlalib.c"
#include "interpolate.c"
#include "matrix.h"


#define SQR(X) X*X


double *GLOBAL_x,*GLOBAL_y,*GLOBAL_xkick1,*GLOBAL_ykick1,*GLOBAL_xkick,*GLOBAL_ykick,*GLOBAL_xkick2,*GLOBAL_ykick2; 
int GLOBAL_m,GLOBAL_n;

/*Definition of the interpolated functions*/
static double Map_x(double x,double y)
{
    double f;
    /*cubic interpolation*/
    /*splin2(GLOBAL_y,GLOBAL_x,GLOBAL_xkick,GLOBAL_xkick2,GLOBAL_n,GLOBAL_m,y,x,&f);*/
    
    /*biliniar interpolation*/
    linint(GLOBAL_y,GLOBAL_x,GLOBAL_xkick,GLOBAL_m,GLOBAL_n,y,x,&f);
    return f;
}

static double Map_y(double x,double y)
{
    double f;
     /*cubic interpolation*/
    /*splin2(GLOBAL_y,GLOBAL_x,GLOBAL_ykick,GLOBAL_ykick2,GLOBAL_m,GLOBAL_n,y,x,&f);*/
        
    /*biliniar interpolation*/
    linint(GLOBAL_y,GLOBAL_x,GLOBAL_ykick,GLOBAL_m,GLOBAL_n,y,x,&f);
    return f;
}

static double Map1_x(double x,double y)
{
    double f;
    /*cubic interpolation*/
    /*splin2(GLOBAL_y,GLOBAL_x,GLOBAL_xkick1,GLOBAL_xkick2,GLOBAL_n,GLOBAL_m,y,x,&f);*/
    
    /*biliniar interpolation*/
    linint(GLOBAL_y,GLOBAL_x,GLOBAL_xkick1,GLOBAL_m,GLOBAL_n,y,x,&f);
    return f;
}

static double Map1_y(double x,double y)
{
    double f;
     /*cubic interpolation*/
    /*splin2(GLOBAL_y,GLOBAL_x,GLOBAL_ykick1,GLOBAL_ykick2,GLOBAL_m,GLOBAL_n,y,x,&f);*/
        
    /*biliniar interpolation*/
    linint(GLOBAL_y,GLOBAL_x,GLOBAL_ykick1,GLOBAL_m,GLOBAL_n,y,x,&f);
    return f;
}


static void markaslost(double *r6)
{	int i;
    
    r6[0] = mxGetNaN();
	
    for(i=1;i<6;i++)
		r6[i] =0;
}


/* Set T1, T2, R1, R2 to NULL pointers to ignore misalignmets*/
void IdKickMapModelPass(double *r, double le, double *xkick1, double *ykick1, double *xkick, double *ykick, double *x, double *y,int n,int m, int Nslice, double *T1, double *T2, double *R1, double *R2, int num_particles)
{	double *r6,f,L1,deltaxp,deltayp,deltaxp1,deltayp1,*limitsptr;   
	int c;
    bool useT1, useT2, useR1, useR2;
    
    /*Act as AperturePass*/
    limitsptr=(double*)mxCalloc(4,sizeof(double));
    limitsptr[0]=x[0];
    limitsptr[1]=x[n-1];
    limitsptr[2]=y[0];
    limitsptr[3]=y[m-1];
    
    
    

     /*globalize*/
    
    /* For cubic interpolation only*/
    
    /*GLOBAL_xkick2=(double*)mxCalloc(n*m,sizeof(double));
	GLOBAL_ykick2=(double*)mxCalloc(n*m,sizeof(double));
    splie2(y,x,xkick,m,n,GLOBAL_xkick2);
    splie2(y,x,ykick,m,n,GLOBAL_ykick2); */

    GLOBAL_x=x;
    GLOBAL_y=y;
    GLOBAL_xkick1=xkick1;
    GLOBAL_ykick1=ykick1;
    GLOBAL_xkick=xkick;
    GLOBAL_ykick=ykick;
    GLOBAL_m=m; /* y used as colums*/
    GLOBAL_n=n; /* x used as rows*/
     
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
	     

	  
    L1=le/(2*Nslice);
	for(c = 0;c<num_particles;c++)
		{	
		    r6 = r+c*6;
            
			if(!mxIsNaN(r6[0]) & mxIsFinite(r6[4]))
			/* 
		       function bend6 internally calculates the square root
			   of the energy deviation of the particle 
			   To protect against DOMAIN and OVERFLOW error, check if the
			   fifth component of the phase spacevector r6[4] is finite
			*/
			{	                
                if(r6[0]<limitsptr[0] || r6[0]>limitsptr[1] || r6[2]<limitsptr[2] || r6[2]>limitsptr[3])
                {
                    markaslost(r6);
                }
                else
                {
                    /* Misalignment at entrance */
                    if(useT1) ATaddvv(r6,T1);
                    if(useR1) ATmultmv(r6,R1);
                    /*Tracking in the main body*/
                    for(m=0; m < Nslice; m++) /* Loop over slices*/			
                    {		  
                        r6 = r+c*6;		  
                        ATdrift6(r6,L1);        
                        if (!mxIsNaN(r6[0])&&!mxIsNaN(r6[2])) 
                        {     /*The kick from IDs varies quadratically, not linearly, with energy.   */
                            deltaxp = (1.0/Nslice)*Map_x(r6[0],r6[2])/(1.0+r6[4]);        
                            deltayp = (1.0/Nslice)*Map_y(r6[0],r6[2])/(1.0+r6[4]); 
                            deltaxp1 = (1.0/Nslice)*Map1_x(r6[0],r6[2]);        
                            deltayp1 = (1.0/Nslice)*Map1_y(r6[0],r6[2]); 
                            r6[1] = r6[1] + deltaxp + deltaxp1; 
                            r6[3] = r6[3] + deltayp + deltayp1;
                        }
                        ATdrift6(r6,L1);	
                    }  
                    /* Misalignment at exit */	
                    if(useR2) ATmultmv(r6,R2);
                    if(useT2) ATaddvv(r6,T2);                    
                }
	        }
		}	
}



/*Just for debugg!! Subtitute routine for a drift*/
static void DriftPass(double *r_in, double le, int num_particles)
/* le - physical length
   r_in - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{	int c, c6;
	double p_norm, NormL;
    mexPrintf("Length: %f\n",le);
    mexPrintf("Initial %e %e %e %e %e %e\n",r_in[0],r_in[1],r_in[2],r_in[3],r_in[4],r_in[5]);
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
    mexPrintf("Final %e %e %e %e %e %e\n",r_in[0],r_in[1],r_in[2],r_in[3],r_in[4],r_in[5]);
}





#ifndef NOMEX
ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
								double *r_in,int num_particles,int mode)

#define NUM_FIELDS_2_REMEMBER 12

{	int *returnptr,n,m;
	int *NewFieldNumbers,fnum,Nslice;
	double le;
	double  *pr1, *pr2, *pt1, *pt2, *xkick, *ykick, *xkick1, *ykick1, *x, *y; 

	
	switch(mode)
		{	case NO_LOCAL_COPY:	/* Get fields by names from MATLAB workspace   */
				{	/* NO_LOCAL_COPY - obsolete since AT1.3 */
				    returnptr = NULL;
				}	break;	

			case MAKE_LOCAL_COPY: 	/* Find field numbers first
									   Save a list of field number in an array
									   and make returnptr point to that array
									*/
				{	
					/* Populate */					
					/*mexPrintf("Make copy\n");*/
                    NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
					
					
					
					fnum = mxGetFieldNumber(ElemData,"Length");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));

					fnum = mxGetFieldNumber(ElemData,"xkick1");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'xkick1' was not found in the element data structure"); 					
					NewFieldNumbers[1] = fnum;
                    xkick1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
                    n = mxGetN(mxGetFieldByNumber(ElemData,0,fnum));
                    m = mxGetM(mxGetFieldByNumber(ElemData,0,fnum));
                    
					fnum = mxGetFieldNumber(ElemData,"ykick1");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'ykick1' was not found in the element data structure"); 					
					NewFieldNumbers[2] = fnum;
                     ykick1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
                     
                    
					fnum = mxGetFieldNumber(ElemData,"xkick");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'xkick' was not found in the element data structure"); 					
					NewFieldNumbers[3] = fnum;
                    xkick = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
                    
					fnum = mxGetFieldNumber(ElemData,"ykick");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'ykick' was not found in the element data structure"); 					
					NewFieldNumbers[4] = fnum;
                     ykick = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
                    
					fnum = mxGetFieldNumber(ElemData,"xtable");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'x' was not found in the element data structure"); 					
					NewFieldNumbers[5] = fnum;
                    x = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
                    
					fnum = mxGetFieldNumber(ElemData,"ytable");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'y' was not found in the element data structure"); 					
					NewFieldNumbers[6] = fnum;
                    y = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
                    
                    
                    
                    					
                    
                    
                    
                    
                    
                    fnum = mxGetFieldNumber(ElemData,"Nslice");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Nslice' was not found in the element data structure"); 					
					NewFieldNumbers[7] = fnum;
                    Nslice = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
                    
					
					/* Optional fields */

				
                    fnum = mxGetFieldNumber(ElemData,"R1");
					NewFieldNumbers[8] = fnum;
					if(fnum<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					

                    fnum = mxGetFieldNumber(ElemData,"R2");
					NewFieldNumbers[9] = fnum;
					if(fnum<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
                    fnum = mxGetFieldNumber(ElemData,"T1");
	                NewFieldNumbers[10] = fnum;
					if(fnum<0)
					    pt1 = NULL;
					else
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
	                fnum = mxGetFieldNumber(ElemData,"T2");
	                NewFieldNumbers[11] = fnum;
					if(fnum<0)
					    pt2 = NULL;
					else
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));

					returnptr = NewFieldNumbers;

				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
									    The second argument ponter to the array of field 
										numbers is previously created with 
										BendLinearPass( ..., MAKE_LOCAL_COPY)
									*/
											
				{	
                    /*mexPrintf("Aqui hi vas?\n");*/
                    le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
                    xkick1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
                    ykick1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
                    xkick = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
                    ykick = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
                    x = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
                    y = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
                    Nslice = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));       
                    n = mxGetN(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
                    m = mxGetM(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
                    
					/* Optional fields */
					if(FieldNumbers[8]<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[8]));
					
					if(FieldNumbers[9]<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[9]));
					
					    
					if(FieldNumbers[10]<0)
					    pt1 = NULL;
					else    
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[10]));
					    
					if(FieldNumbers[11]<0)
					    pt2 = NULL;
					else 
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[11]));



					returnptr = FieldNumbers;

				}	break;

		}      
    
    
    /*DriftPass(r_in,le,num_particles);*/
    IdKickMapModelPass(r_in, le,xkick1,ykick1,xkick,ykick,x,y,n,m,Nslice,  
							pt1, pt2, pr1, pr2, num_particles);
	
    return(returnptr);
}







void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	int m,n,m_map,n_map,Nslice;
	double *r_in,le,*xkick,*ykick,*xkick1,*ykick1,*x,*y; 
	double  *pr1, *pr2, *pt1, *pt2;
	mxArray *tmpmxptr;
  
	if(nrhs)
	{
        /* ALLOCATE memory for the output array of the same size as the input */
	    m = mxGetM(prhs[1]);
	    n = mxGetN(prhs[1]);
	    if(m!=6) 
	    {mexErrMsgTxt("Second argument must be a 6 x N matrix");}	
	    
	    tmpmxptr = mxGetField(prhs[0],0,"Length");
	    if(tmpmxptr)
	        le = mxGetScalar(tmpmxptr);
	    else
	        mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
	    
        tmpmxptr = mxGetField(prhs[0],0,"xkick1");
	    if(tmpmxptr)
        {
            xkick1 = mxGetPr(tmpmxptr);
            n_map = mxGetN(tmpmxptr);        
            m_map = mxGetM(tmpmxptr);
        }
	    else
	        mexErrMsgTxt("Required field 'xkick1' was not found in the element data structure");    	    
        
 	    tmpmxptr = mxGetField(prhs[0],0,"ykick1");
	    if(tmpmxptr)
	        ykick1 = mxGetPr(tmpmxptr);
	    else
	        mexErrMsgTxt("Required field 'ykick1' was not found in the element data structure");   
        
        
	    tmpmxptr = mxGetField(prhs[0],0,"xkick");
	    if(tmpmxptr)
        {
            xkick = mxGetPr(tmpmxptr);
        }
	    else
	        mexErrMsgTxt("Required field 'xkick' was not found in the element data structure");    	    
        
 	    tmpmxptr = mxGetField(prhs[0],0,"ykick");
	    if(tmpmxptr)
	        ykick = mxGetPr(tmpmxptr);
	    else
	        mexErrMsgTxt("Required field 'ykick' was not found in the element data structure");           
        
 	    tmpmxptr = mxGetField(prhs[0],0,"xtable");
	    if(tmpmxptr)
	        x = mxGetPr(tmpmxptr);
	    else
	        mexErrMsgTxt("Required field 'xtable' was not found in the element data structure");           
        
	    tmpmxptr = mxGetField(prhs[0],0,"ytable");
	    if(tmpmxptr)
	        y = mxGetPr(tmpmxptr);
	    else
	        mexErrMsgTxt("Required field 'ytable' was not found in the element data structure");    
        
  	    tmpmxptr = mxGetField(prhs[0],0,"Nslice");
	    if(tmpmxptr)
	        Nslice = (int)mxGetScalar(tmpmxptr);
	    else
	        mexErrMsgTxt("Required field 'Nslice' was not found in the element data structure");       
        
        
        
        /*Optional fields*/ 
        
        
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

		IdKickMapModelPass(r_in, le, xkick1,ykick1, xkick,ykick,x,y,n_map,m_map,Nslice,  
							pt1, pt2, pr1, pr2, n);
	}
	else                             
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(8,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("xkick1"));
	    mxSetCell(plhs[0],2,mxCreateString("ykick1"));
	    mxSetCell(plhs[0],3,mxCreateString("xkick"));
	    mxSetCell(plhs[0],4,mxCreateString("ykick"));
	    mxSetCell(plhs[0],5,mxCreateString("xtable"));
        mxSetCell(plhs[0],6,mxCreateString("ytable"));
        mxSetCell(plhs[0],7,mxCreateString("Nslice"));
                        
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(4,1);
	        mxSetCell(plhs[1],0,mxCreateString("T1"));
	        mxSetCell(plhs[1],1,mxCreateString("T2"));
	        mxSetCell(plhs[1],2,mxCreateString("R1"));
	        mxSetCell(plhs[1],3,mxCreateString("R2"));
	    }
	}
 }
#endif
