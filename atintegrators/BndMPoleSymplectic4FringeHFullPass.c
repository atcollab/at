#include "mex.h"
#include<math.h>
#include "elempass.h"
#include "atlalib.c"
#include "atphyslib.c"

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656


#define SQR(X) ((X)*(X))



void bndthinkick(double* r, double* A, double* B, double L, double irho, int max_order)

/***************************************************************************** 
Calculate multipole kick in a curved elemrnt (bending magnet)
The reference coordinate system  has the curvature given by the inverse 
(design) radius irho.
IMPORTANT !!!
The magnetic field Bo that provides this curvature MUST NOT be included in the dipole term
PolynomB[1](MATLAB notation)(C: B[0] in this function) of the By field expansion

The kick is given by

           e L      L delta      L x
theta  = - --- B  + -------  -  -----  , 
     x     p    y     rho           2
            0                    rho

         e L
theta  = --- B
     y    p   x
           0


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
   
   /* recursively calculate the local transvrese magnetic field
    * Bx = ReSum, By = ImSum
    */
   for(i=max_order-1;i>0;i--)
   {	ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
        ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
        ReSum = ReSumTemp;
   } 
   ReSumTemp = ReSum*r[0] - ImSum*r[2];
   ImSum = ImSum*r[0] +  ReSum*r[2] + A[0];
   ReSum = ReSumTemp;
   
   r[1] -=  L*ReSum*(1+r[0]*irho);
   r[3] +=  L*ImSum*(1+r[0]*irho);  
   
}

void BndMPoleSymplectic4Pass(double *r, double le, double irho, double *A, double *B,
					int max_order, int num_int_steps,
					double entrance_angle, 	double exit_angle,
					double fint1, double fint2, double gap,
					double *T1, double *T2,	
					double *R1, double *R2, int num_particles)


{	int c,m;	
	double *r6;   
	double SL, L1, L2, K1, K2;
	bool useT1, useT2, useR1, useR2, useFringe1, useFringe2;
	
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
	    
	/* if either is 0 - do not calculate fringe effects */    
    if( fint1==0 || gap==0) 
	    useFringe1 = false;
	else 
	    useFringe1=true;  
	
	if( fint2==0 || gap==0) 
	    useFringe2 = false;
	else 
	    useFringe2=true;  
	
	    
	for(c = 0;c<num_particles;c++)	/* Loop over particles  */
			{	r6 = r+c*6;	
			    if(!mxIsNaN(r6[0]))
			    {
					

                    
					/*  misalignment at entrance  */
					if(useT1)
			            ATaddvv(r6,T1);
			        if(useR1)
			            ATmultmv(r6,R1);    
                    
                    /* edge focus */

                    Rotation_y(r6, entrance_angle);
                    if(useFringe1){                     
                        QuadFringe_step(r6, B[1], 1.0,DRIFT1);
                        BendFringe_FINT(r6, (irho+B[0]), 1.0,fint1,gap,KICK1);
                        QuadFringe_step(r6, B[1], 1.0,DRIFT2);
                        BendFringe_FINT(r6, (irho+B[0]), 1.0,fint1,gap,KICK2);
                        QuadFringe_step(r6, B[1], 1.0,DRIFT2);
                        BendFringe_FINT(r6, (irho+B[0]), 1.0,fint1,gap,KICK1);
                        QuadFringe_step(r6, B[1], 1.0,DRIFT1);
                    }
                    else {
                        QuadFringe_step(r6, B[1], 1.0,DRIFT1);
                        BendFringe(r6, (irho+B[0]), 1.0,KICK1);
                        QuadFringe_step(r6, B[1], 1.0,DRIFT2);
                        BendFringe(r6, (irho+B[0]), 1.0,KICK2);
                        QuadFringe_step(r6, B[1], 1.0,DRIFT2);
                        BendFringe(r6, (irho+B[0]), 1.0,KICK1);
                        QuadFringe_step(r6, B[1], 1.0,DRIFT1);
                    }
                    Ideal_Wedge(r6, -entrance_angle, irho);

					/* integrator */
					for(m=0; m < num_int_steps; m++) /* Loop over slices*/			
						{		r6 = r+c*6;	
                                /*printf("At bend s0 x %1.8f px %1.8f y %1.8f py %1.8f d %1.8f l %1.8f\n",r6[0],r6[1],r6[2],r6[3],r6[4],r6[5]);*/
								AT_H_Full_Bend_y_Drift(r6,L1,irho);
           					    bndthinkick(r6, A, B, K1,irho, max_order);
								AT_H_Full_Bend_y_Drift(r6,L2,irho);
           					    bndthinkick(r6, A, B, K2,irho, max_order);
								AT_H_Full_Bend_y_Drift(r6,L2,irho);
		     					bndthinkick(r6, A, B,  K1,irho, max_order);
								AT_H_Full_Bend_y_Drift(r6,L1,irho);	
						}  
                    /* edge focus */

                    Ideal_Wedge(r6, -exit_angle, irho);
                    if(useFringe2){                     
                        QuadFringe_step(r6, B[1], -1.0,DRIFT1);
                        BendFringe_FINT(r6, (irho+B[0]), -1.0,fint2,gap,KICK1);
                        QuadFringe_step(r6, B[1], -1.0,DRIFT2);
                        BendFringe_FINT(r6, (irho+B[0]), -1.0,fint2,gap,KICK2);
                        QuadFringe_step(r6, B[1], -1.0,DRIFT2);
                        BendFringe_FINT(r6, (irho+B[0]), -1.0,fint2,gap,KICK1);
                        QuadFringe_step(r6, B[1], -1.0,DRIFT1);
                    }
                    else {
                        QuadFringe_step(r6, B[1], -1.0,DRIFT1);
                        BendFringe(r6, (irho+B[0]), -1.0,KICK1);
                        QuadFringe_step(r6, B[1], -1.0,DRIFT2);
                        BendFringe(r6, (irho+B[0]), -1.0,KICK2);
                        QuadFringe_step(r6, B[1], -1.0,DRIFT2);
                        BendFringe(r6, (irho+B[0]), -1.0,KICK1);
                        QuadFringe_step(r6, B[1], -1.0,DRIFT1);
                    }
                    Rotation_y(r6, exit_angle);

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

#define NUM_FIELDS_2_REMEMBER 15


{	double *A , *B;
	double  *pr1, *pr2, *pt1, *pt2, fint1, fint2, gap;   
	double entrance_angle, exit_angle;

	int max_order, num_int_steps;
	double le,ba,irho;
	int *returnptr;
	int *NewFieldNumbers, fnum;

	
	switch(mode)
		{   case MAKE_LOCAL_COPY: 	/* Find field numbers first
										Save a list of field number in an array
										and make returnptr point to that array
									*/
				{	
					/* Allocate memory for integer array of 
					    field numbers for faster future reference
					*/
		
					NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));

					/* Populate */
					
					
					
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
					
					
					fnum = mxGetFieldNumber(ElemData,"BendingAngle");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'BendingAngle' was not found in the element data structure"); 
					NewFieldNumbers[5] = fnum;
					ba = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					
					
					
	                fnum = mxGetFieldNumber(ElemData,"EntranceAngle");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'EntranceAngle' was not found in the element data structure"); 
					NewFieldNumbers[6] = fnum;
					entrance_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
	                
	                fnum = mxGetFieldNumber(ElemData,"ExitAngle");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'ExitAngle' was not found in the element data structure"); 
					NewFieldNumbers[7] = fnum;
					exit_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					
					fnum = mxGetFieldNumber(ElemData,"FringeInt1");/* Optional field FringeInt */
                    NewFieldNumbers[8] = fnum;
					if(fnum<0) 
					    fint1 = 0;
					else
					    fint1 = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					    
					    
					fnum = mxGetFieldNumber(ElemData,"FringeInt2");/* Optional field FringeInt */
                    NewFieldNumbers[9] = fnum;
					if(fnum<0) 
					    fint2 = 0;
					else
					    fint2 = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					fnum = mxGetFieldNumber(ElemData,"FullGap");
					NewFieldNumbers[10] = fnum;
					if(fnum<0) 
					    gap = 0;
					else
					    gap = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
				
                    fnum = mxGetFieldNumber(ElemData,"R1");
					NewFieldNumbers[11] = fnum;
					if(fnum<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					

					fnum = mxGetFieldNumber(ElemData,"R2");
					NewFieldNumbers[12] = fnum;
					if(fnum<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
                    fnum = mxGetFieldNumber(ElemData,"T1");
	                NewFieldNumbers[13] = fnum;
					if(fnum<0)
					    pt1 = NULL;
					else
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
	                
	                fnum = mxGetFieldNumber(ElemData,"T2");
	                NewFieldNumbers[14] = fnum;
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
					ba = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
					entrance_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
					exit_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));
					
					/* Optional fields */
					
					if(FieldNumbers[8]<0) 
					    fint1 = 0;
					else
					    fint1 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[8]));
					
					    
					if(FieldNumbers[9]<0) 
					    fint2 = 0;
					else
					    fint2 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[9]));
					    
					if(FieldNumbers[10]<0) 
					    gap = 0;
					else
					gap = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[10]));
					
					/* Optional fields */
					if(FieldNumbers[11]<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[11]));
					
					if(FieldNumbers[12]<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[12]));
					
					    
					if(FieldNumbers[13]<0)
					    pt1 = NULL;
					else    
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[13]));
					    
					if(FieldNumbers[14]<0)
					    pt2 = NULL;
					else 
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[14]));
					
			
					
					returnptr = FieldNumbers;
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function BndMPoleSymplectic4Pass\n");
				}
		}


	irho = ba/le;
	
	BndMPoleSymplectic4Pass(r_in, le, irho, A, B, max_order, num_int_steps, 
								entrance_angle, exit_angle, fint1, fint2, gap, pt1, pt2, pr1, pr2, num_particles);
	

	
	return(returnptr);

}


 






void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	int m,n;
	double *r_in;
	double le, ba, *A, *B;  
	double irho;
	int max_order, num_int_steps;
	double entrance_angle, exit_angle ;
	double  *pr1, *pr2, *pt1, *pt2, fint1, fint2, gap;  
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
		
		
    irho = ba/le;
    plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	BndMPoleSymplectic4Pass(r_in, le, irho, A, B, max_order, num_int_steps, 
								entrance_angle, exit_angle, fint1, fint2, gap, pt1, pt2, pr1, pr2, n);

	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(8,1);
	    
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("BendingAngle"));
	    mxSetCell(plhs[0],2,mxCreateString("EntranceAngle"));
	    mxSetCell(plhs[0],3,mxCreateString("ExitAngle"));
        mxSetCell(plhs[0],4,mxCreateString("PolynomA"));
	    mxSetCell(plhs[0],5,mxCreateString("PolynomB"));
	    mxSetCell(plhs[0],6,mxCreateString("MaxOrder"));
	    mxSetCell(plhs[0],7,mxCreateString("NumIntSteps"));	 	    
	    
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(7,1);
	        mxSetCell(plhs[1],0,mxCreateString("FullGap"));
	        mxSetCell(plhs[1],1,mxCreateString("FringeInt1"));
	        mxSetCell(plhs[1],2,mxCreateString("FringeInt2"));
	        mxSetCell(plhs[1],3,mxCreateString("T1"));
	        mxSetCell(plhs[1],4,mxCreateString("T2"));
	        mxSetCell(plhs[1],5,mxCreateString("R1"));
	        mxSetCell(plhs[1],6,mxCreateString("R2"));
	    }
	}



}



