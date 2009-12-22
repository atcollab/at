/* BendLinearPass.c
   Accelerator Toolbox 
   Revision 6/26/00
   A.Terebilo terebilo@ssrl.slac.stanford.edu
*/

#include "mex.h"
#include "elempass.h"
#include "atlalib.c"
#include "matrix.h"
#include <math.h>

#define SQR(X) X*X

void edge(double* r, double inv_rho, double edge_angle)
{	double psi = inv_rho*tan(edge_angle);
	r[1]+=r[0]*psi;
	r[3]-=r[2]*psi;
}


void edge_fringe(double* r, double inv_rho, double edge_angle, double fint, double gap)
{ double fx = inv_rho*tan(edge_angle);
  double psi_bar = edge_angle-inv_rho*gap*fint*(1+sin(edge_angle)*sin(edge_angle))/cos(edge_angle);
	double fy = inv_rho*tan(psi_bar);
	r[1]+=r[0]*fx;
	r[3]-=r[2]*fy;
}



void bend6(double* r, double L, double b_angle, double grd, double ByError)
{	
	double M12,M21,M34,M43,MVD,MHD;  /*  non-0 elements of transfer matrix */
	double x, xpr,y,ypr,delta;
	
	
	
	double sqrtG1, sqrtG2, arg1, arg2;
	double Kx = b_angle/L; /* Curvature of the design trajectory */
	double p_norm = 1/(1+r[4]);
	double G1 = (Kx*Kx+grd)*p_norm;
	double G2 = -grd*p_norm;

	
	/* Horizontal Transverse Matrix */
	if(G1==0)		/* Special case Kx^2 + grd = 0 */
 
		{	MHD = 1;				
			M12 = L;
			M21 = 0;	
		}
	else		
		{	if(G1 > 0)
				{	sqrtG1 = sqrt(G1);
					arg1 = L*sqrtG1;
					MHD = cos(arg1);				
					M12 = sin(arg1)/sqrtG1;
					M21 = -sin(arg1)*sqrtG1;	
				}
			else
				{	sqrtG1 = sqrt(-G1);
					arg1 = L*sqrtG1;
					MHD = cosh(arg1); 				
					M12 = sinh(arg1)/sqrtG1;
					M21 = sinh(arg1)*sqrtG1;	
				}
		}

	

	/*  Vertical Transverse Matrix */
	
	if(G2==0) /*  No gradient - vertical motion is a drift  */
	
		{	MVD = 1;				
			M34 = L;
			M43 = 0;	
		}
	else	
		{	if(G2 > 0)	/* Vertical focusing */
				{	sqrtG2 = sqrt(G2);
					arg2 = L*sqrtG2;
					MVD = cos(arg2);;				
					M34 = sin(arg2)/sqrtG2;
					M43 = -sin(arg2)*sqrtG2;	
				}
			else		/*  Vertical defocusing	*/
				{	sqrtG2 = sqrt(-G2);
					arg2 = L*sqrtG2;
					MVD = cosh(arg2); 				
					M34 = sinh(arg2)/sqrtG2;
					M43 = sinh(arg2)*sqrtG2;	
				}
		}

	x   = r[0];
	xpr = r[1]*p_norm;
	y   = r[2];
	ypr = r[3]*p_norm;
	delta = r[4]; 

	r[0]=  MHD*x + M12*xpr ;
	r[1]= (M21*x + MHD*xpr)/p_norm; 
	
	if(G1==0)	
		{	r[0]+= (delta*p_norm-ByError)*L*L*Kx/2;
			r[1]+= (delta*p_norm-ByError)*L*Kx/p_norm ;
		}
	else
		{	if(G1>0)
				{	r[0]+= (delta*p_norm-ByError)*(1-cos(arg1))*Kx/G1;
					r[1]+= (delta*p_norm-ByError)*sin(arg1)*Kx/(sqrtG1*p_norm) ;
				}	
			else
				{	r[0]+= (delta*p_norm-ByError)*(1-cosh(arg1))*Kx/G1;
					r[1]+= (delta*p_norm-ByError)*sinh(arg1)*Kx/(sqrtG1*p_norm) ;
				}
		}

	r[2]=  MVD*y + M34*ypr;
	r[3]= (M43*y + MVD*ypr)/p_norm ;
		
	r[5]+= xpr*xpr*(L+MHD*M12)/4;
   	
	if (G1==0) {
	    /* Do nothing */   
	} else {
    	r[5]+= (L-MHD*M12)*(x*x*G1+(delta*p_norm-ByError)*(delta*p_norm-ByError)*Kx*Kx/G1 
										-2*x*Kx*(delta*p_norm-ByError))/4;

    	r[5]+= M12*M21*( x*xpr - xpr*(delta*p_norm-ByError)*Kx/G1)/2;
    	
    	r[5]+= Kx*x*M12  +   xpr*(1-MHD)*Kx/G1   +   (delta*p_norm-ByError)*(L-M12)*Kx*Kx/G1;
    }



	r[5]+= ((L-MVD*M34)*y*y*G2 + ypr*ypr*(L+MVD*M34))/4;
   	
	r[5]+= M34*M43*x*xpr/2;	



}



void BendLinearPass(double *r, double le, double grd ,double ba, double bye,	
					    double entrance_angle, double exit_angle,
						double fint1, double fint2, double gap,
						double *T1, double *T2, double *R1, double *R2,
						int num_particles)
/* Set T1, T2, R1, R2 to NULL pointers to ignore misalignmets
   Set fint OR gap are 0 to ignore fringe effects
   Set bye to 0 to ignore ByError*/

						
{  double *r6;   
	int c;
    bool useT1, useT2, useR1, useR2, useFringe1, useFringe2, useByError;
	
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
	    
	/* if either is 0 - do not calculate fringe effects */    
    if( fint1==0 || gap==0) 
	    useFringe1 = false;
	else 
	    useFringe1=true;  
	
	if( fint2==0 || gap==0) 
	    useFringe2 = false;
	else 
	    useFringe2=true;  
	    
	if(bye==0) 
	    useByError = false;
	else 
	    useByError=true;  

	    
	
	for(c = 0;c<num_particles;c++)
		{	
		    r6 = r+c*6;

			if(!mxIsNaN(r6[0]) || !mxIsFinite(r6[4]))
			/* 
		       function bend6 internally calculates the square root
			   of the energy deviation of the particle 
			   To protect against DOMAIN and OVERFLOW error, check if the
			   fifth component of the phase spacevector r6[4] is finite
			*/
			{	/* Misalignment at entrance */
	            if(useT1)
			        ATaddvv(r6,T1);
			    if(useR1)
			        ATmultmv(r6,R1);
			    
			    if(useFringe1)
			        edge_fringe(r6, ba/le, entrance_angle,fint1,gap);
			    else
			        edge(r6, ba/le, entrance_angle);		 
				  
			    bend6(r6, le, ba, grd, bye);
			    
			    if(useFringe2)
			        edge_fringe(r6, ba/le, exit_angle,fint2,gap);
			    else
			        edge(r6, ba/le, exit_angle);
    
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

#define NUM_FIELDS_2_REMEMBER 13

{	int *returnptr;
	int *NewFieldNumbers,fnum;
	double le, bye, ba, grd, entrance_angle, exit_angle, fint1, fint2 , gap ;
	double  *pr1, *pr2, *pt1, *pt2; 

	
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
					NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
					
					
					
					fnum = mxGetFieldNumber(ElemData,"Length");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));
					
                    
					
					
					
					fnum = mxGetFieldNumber(ElemData,"BendingAngle");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'BendingAngle' was not found in the element data structure"); 					
					NewFieldNumbers[1] = fnum;
                    ba = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
                    fnum = mxGetFieldNumber(ElemData,"EntranceAngle");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'EntranceAngle' was not found in the element data structure"); 					
					NewFieldNumbers[2] = fnum;
					entrance_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));

					
					fnum = mxGetFieldNumber(ElemData,"ExitAngle");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'ExitAngle' was not found in the element data structure"); 					
					NewFieldNumbers[3] = fnum;
					exit_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					/* Optional fields */
					
					fnum = mxGetFieldNumber(ElemData,"K"); /* Optional field K */
					NewFieldNumbers[4] = fnum;
					if(fnum<0) 
					    grd = 0;
					else
					    grd = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					fnum = mxGetFieldNumber(ElemData,"ByError"); /* Optional field ByError */
                    NewFieldNumbers[5] = fnum;
					if(fnum<0) 
					    bye = 0;
					else
					    bye = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));    
					
					fnum = mxGetFieldNumber(ElemData,"FringeInt1");/* Optional field FringeInt */
                    NewFieldNumbers[6] = fnum;
					if(fnum<0) 
					    fint1 = 0;
					else
					    fint1 = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					    
					fnum = mxGetFieldNumber(ElemData,"FringeInt2");/* Optional field FringeInt */
                    NewFieldNumbers[7] = fnum;
					if(fnum<0) 
					    fint2 = 0;
					else
					    fint2 = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					fnum = mxGetFieldNumber(ElemData,"FullGap");
					NewFieldNumbers[8] = fnum;
					if(fnum<0) 
					    gap = 0;
					else
					    gap = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
				
                    fnum = mxGetFieldNumber(ElemData,"R1");
					NewFieldNumbers[9] = fnum;
					if(fnum<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					

					fnum = mxGetFieldNumber(ElemData,"R2");
					NewFieldNumbers[10] = fnum;
					if(fnum<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
                    fnum = mxGetFieldNumber(ElemData,"T1");
	                NewFieldNumbers[11] = fnum;
					if(fnum<0)
					    pt1 = NULL;
					else
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
	                
	                fnum = mxGetFieldNumber(ElemData,"T2");
	                NewFieldNumbers[12] = fnum;
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
				    le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					ba = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
					entrance_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
					exit_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
					
					
					/* Optional fields */
					if(FieldNumbers[4]<0) 
					    grd = 0;
					else
					    grd = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
					
					if(FieldNumbers[5]<0) 
					    bye = 0;
					else
					    bye = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
					    
					if(FieldNumbers[6]<0) 
					    fint1 = 0;
					else
					    fint1 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
					    
					if(FieldNumbers[7]<0) 
					    fint2 = 0;
					else
					    fint2 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));
					
					if(FieldNumbers[8]<0) 
					    gap = 0;
					else
					gap = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[8]));
					
					/* Optional fields */
					if(FieldNumbers[9]<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[9]));
					
					if(FieldNumbers[10]<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[10]));
					
					    
					if(FieldNumbers[11]<0)
					    pt1 = NULL;
					else    
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[11]));
					    
					if(FieldNumbers[12]<0)
					    pt2 = NULL;
					else 
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[12]));



					returnptr = FieldNumbers;

				}	break;

		}
		
		BendLinearPass(r_in, le, grd , ba, bye,	
							entrance_angle, exit_angle, fint1, fint2, gap,  
							pt1, pt2, pr1, pr2, num_particles);


		return(returnptr);
}







void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	int m,n;
	double *r_in, le, grd, ba, bye,	entrance_angle, exit_angle, fint1, fint2 , gap; 
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
	    
	    tmpmxptr = mxGetField(prhs[0],0,"ByError");
	    if(tmpmxptr)
	        bye = mxGetScalar(tmpmxptr);
	    else
	        bye = 0; 
	    
	    tmpmxptr = mxGetField(prhs[0],0,"BendingAngle");
	    if(tmpmxptr)
	        ba = mxGetScalar(tmpmxptr);
	    else
	        mexErrMsgTxt("Required field 'BendingAngle' was not found in the element data structure"); 
	    
	    tmpmxptr = mxGetField(prhs[0],0,"EntranceAngle");
	    if(tmpmxptr)
	        entrance_angle = mxGetScalar(tmpmxptr);
	    else
	        mexErrMsgTxt("Required field 'EntranceAngle' was not found in the element data structure"); 
	    
	    
	    tmpmxptr = mxGetField(prhs[0],0,"ExitAngle");
	    if(tmpmxptr)
	        exit_angle = mxGetScalar(tmpmxptr);
	    else
	        mexErrMsgTxt("Required field 'ExitAngle' was not found in the element data structure");     
	    
	    
	    tmpmxptr = mxGetField(prhs[0],0,"K");
	    if(tmpmxptr)
	        grd = mxGetScalar(tmpmxptr);
	    else
	        grd = 0;	    
	    
	    
	    tmpmxptr = mxGetField(prhs[0],0,"FringeInt1");
	    if(tmpmxptr)
	        fint1 = mxGetScalar(tmpmxptr);
	    else
	        fint1 = 0;
	        
	        
	    tmpmxptr = mxGetField(prhs[0],0,"FringeInt2");
	    if(tmpmxptr)
	        fint2 = mxGetScalar(tmpmxptr);
	    else
	        fint2 = 0;
	    
	    tmpmxptr = mxGetField(prhs[0],0,"FullGap");
	    if(tmpmxptr)
	        gap = mxGetScalar(tmpmxptr);
	    else
	        gap = 0;
	    
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

 	    BendLinearPass(r_in, le, grd, ba, bye, entrance_angle, exit_angle,fint1, fint2, gap, 
												pt1, pt2, pr1, pr2, n);
	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(4,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("BendingAngle"));
	    mxSetCell(plhs[0],2,mxCreateString("EntranceAngle"));
	    mxSetCell(plhs[0],3,mxCreateString("ExitAngle"));
	    
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(9,1);
	        mxSetCell(plhs[1],0,mxCreateString("K"));
	        mxSetCell(plhs[1],1,mxCreateString("ByError"));
	        mxSetCell(plhs[1],2,mxCreateString("FullGap"));
	        mxSetCell(plhs[1],3,mxCreateString("FringeInt1"));
	        mxSetCell(plhs[1],4,mxCreateString("FringeInt2"));
	        mxSetCell(plhs[1],5,mxCreateString("T1"));
	        mxSetCell(plhs[1],6,mxCreateString("T2"));
	        mxSetCell(plhs[1],7,mxCreateString("R1"));
	        mxSetCell(plhs[1],8,mxCreateString("R2"));
	    }
	}




 }
