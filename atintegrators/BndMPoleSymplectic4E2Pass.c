#include "at.h"
#include "atlalib.c"
#include "atphyslib.c"

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

/*
 This code was modified from the original BndMPoleSymplectic4Pass.c of AT to correctly integrate the Hamiltonian in 
 the curvilinear coordinate system of the dipole and to include the second order Transport map of the fringe field. 
 New version created by Xiaobiao Huang in March 2009, in final verified version in August 2009.

 */
#define SQR(X) ((X)*(X))

void edge_fringe2A(double* r, double inv_rho, double edge_angle, double fint, double gap,double h1,double K1);
void edge_fringe2B(double* r, double inv_rho, double edge_angle, double fint, double gap,double h2,double K1);
void ATmultmv(double *r, const double* A);
void ATaddvv(double *r, const double *dr);

/*original kick function by Andrei Terebilo*/
static void bndthinkick0(double* r, double* A, double* B, double L, double irho, int max_order)

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
	   Bx = ReSum, By = ImSum
	*/
	for(i=max_order-1;i>=0;i--)
		{	ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
			ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
			ReSum = ReSumTemp;
		}
	
  	
	r[1] -=  L*(ReSum-(r[4]-r[0]*irho)*irho);
	r[3] +=  L*ImSum;
	r[5] +=  L*irho*r[0]; /* pathlength */


}

static void bndthinkick(double* r, double* A, double* B, double L, double h, int max_order)
/*****************************************************************************
(1) PolynomA is neglected.
(2) The vector potential is expanded up to 4th order of x and y. 
(3) Coefficients in PolynomB higher than 4th order is treated as if they are on straight geometry.
(4) The Hamiltonian is H2 = - h x delta - (1+h x)As/Brho-B0 x/Brho      
*/
{
    int i;
	double ReSum = 0; /*B[max_order];*/
 	double ImSum = 0; /*A[max_order];*/

	double ReSumTemp;
	double K1,K2;
 
	K1 = B[1];
	K2 = (max_order>=2) ? B[2] : 0;

    ReSum = B[max_order];
    for(i=max_order-1;i>=0;i--) {
        ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
        ImSum = ImSum*r[0] +  ReSum*r[2] ;
        ReSum = ReSumTemp;
    }

    r[1] -=  L*(-h*r[4] + ReSum + h*(h*r[0]+K1*(r[0]*r[0]-0.5*r[2]*r[2])+K2*(r[0]*r[0]*r[0]-4.0/3.0*r[0]*r[2]*r[2]))    );
    r[3] +=  L*(ImSum+h*(K1*r[0]*r[2]+4.0/3.0*K2*r[0]*r[0]*r[2]+(h/6.0*K1-K2/3.0)*r[2]*r[2]*r[2])) ;
    r[5] +=  L*h*r[0]; /* pathlength */

}



/* the pseudo-drift element described by Hamiltonian H1 = (1+hx) (px^2+py^2)/2(1+delta),     */
static void ATbendhxdrift6(double* r, double L,double h)
{
	double hs = h*L;
	double i1pd = 1.0/(1+r[4]);
	double x=r[0],px=r[1],py=r[3];

	r[0] += (1+h*x)*px*i1pd*L+1/4.*hs*L*(px*px-py*py)*i1pd*i1pd; /* (1.0/h+x)*((1.0+hs*px*i1pd/2.)*(1.0+hs*px*i1pd/2.)-(hs*py*i1pd/2.)*(hs*py*i1pd/2.))-1./h;*/
	r[1] -= hs*(px*px+py*py)*i1pd/2.0;
	
	r[2]+= (1.0+h*x)*i1pd*py*L*(1.+px*hs/2.0);
	r[5]+= (1.0+h*x)*i1pd*i1pd*L/2.0*(px*px+py*py);

}

void BndMPoleSymplectic4E2Pass(double *r, double le, double irho, double *A, double *B,
					int max_order, int num_int_steps,
					double entrance_angle, 	double exit_angle,
					double fint1, double fint2, double gap,double h1,double h2,
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
				 	if(useFringe1)
			    {
			       edge_fringe2A(r6, irho, entrance_angle,fint1,gap,h1,B[1]);
				   /* edge_fringe(r6, irho, entrance_angle,fint1,gap);*/
			    }
			    else
			    {
			       /* edge(r6, irho, entrance_angle); */
			       edge_fringe2A(r6, irho, entrance_angle,0,0,h1,B[1]);

			    }
			          
				 	
					/* integrator */
					for(m=0; m < num_int_steps; m++) /* Loop over slices*/			
						{		r6 = r+c*6;	

								ATbendhxdrift6(r6,L1,irho);
           					    bndthinkick(r6, A, B, K1, irho, max_order);

								ATbendhxdrift6(r6,L2,irho);
           					    bndthinkick(r6, A, B, K2, irho, max_order);
								ATbendhxdrift6(r6,L2,irho);

								bndthinkick(r6, A, B,  K1, irho, max_order);
								ATbendhxdrift6(r6,L1,irho);
						}  
					
					/* edge focus */
				 	if(useFringe2)
			    {
			        edge_fringe2B(r6, irho, exit_angle,fint2,gap,h2,B[1]);
				/*	edge_fringe(r6, irho, exit_angle,fint2,gap);  */
			    }
			    else
			    {
			       /*  edge(r6, irho, exit_angle);	*/
			       edge_fringe2B(r6, irho, exit_angle,0,0,h2,B[1]);
			    }


					 /* Misalignment at exit */	
			        if(useR2)
			            ATmultmv(r6,R2);
		            if(useT2)   
			            ATaddvv(r6,T2);

					
				}


			}
}

#ifdef PYAT

#include "pyutils.c"

int atpyPass(double *rin, int num_particles, PyObject *element, struct parameters *param)
{
    double length = py_get_double(element, "Length", false);
    double bending_angle = py_get_double(element, "BendingAngle", false);
    double entrance_angle = py_get_double(element, "EntranceAngle", false);
    double exit_angle = py_get_double(element, "ExitAngle", false);
    double gap = py_get_double(element, "FullGap", true);
    double fint1 = py_get_double(element, "FringeInt1", true);
    double fint2 = py_get_double(element, "FringeInt2", true);
    long max_order = py_get_long(element, "MaxOrder", true);
    long num_int_steps = py_get_long(element, "NumIntSteps", true);
    double *t1 = numpy_get_double_array(element, "T1", true);
    double *t2 = numpy_get_double_array(element, "T2", true);
    double *r1 = numpy_get_double_array(element, "R1", true);
    double *r2 = numpy_get_double_array(element, "R2", true);
    double *polyA = numpy_get_double_array(element, "PolynomA", true);
    double *polyB = numpy_get_double_array(element, "PolynomB", true);
    double irho = bending_angle / length;
    if (PyErr_Occurred()) {
        return -1;
    } else {
        BndMPoleSymplectic4E2Pass(rin, length, irho, polyA, polyB, max_order, num_int_steps, entrance_angle, exit_angle, fint1, fint2, gap, 0, 0, t1, t2, r1, r2, num_particles);
    return 0;
    }
}

#endif /*PYAT*/

#ifdef MATLAB_MEX_FILE

#include "elempass.h"

ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
								double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 17


{	double *A , *B;
	double  *pr1, *pr2, *pt1, *pt2, fint1, fint2, gap,h1,h2;   
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
					  /*  mexErrMsgTxt("Required field 'PolynomA' was not found in the element data structure");*/
					  A =NULL;
					else
					{
						NewFieldNumbers[0] = fnum;
						A = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					}
					
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
					
					fnum = mxGetFieldNumber(ElemData,"H1");
					NewFieldNumbers[11] = fnum;
					if(fnum<0) 
					    h1 = 0;
					else
					    h1 = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));

					fnum = mxGetFieldNumber(ElemData,"H2");
					NewFieldNumbers[12] = fnum;
					if(fnum<0) 
					    h2 = 0;
					else
					    h2 = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));

				
                    fnum = mxGetFieldNumber(ElemData,"R1");
					NewFieldNumbers[13] = fnum;
					if(fnum<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					

					fnum = mxGetFieldNumber(ElemData,"R2");
					NewFieldNumbers[14] = fnum;
					if(fnum<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
                    fnum = mxGetFieldNumber(ElemData,"T1");
	                NewFieldNumbers[15] = fnum;
					if(fnum<0)
					    pt1 = NULL;
					else
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
	                
	                fnum = mxGetFieldNumber(ElemData,"T2");
	                NewFieldNumbers[16] = fnum;
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

					if(FieldNumbers[11]<0) 
					    h1 = 0;
					else
					h1 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[11]));

					if(FieldNumbers[12]<0) 
					    h2 = 0;
					else
					h2 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[12]));
					
					/* Optional fields */
					if(FieldNumbers[13]<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[13]));
					
					if(FieldNumbers[14]<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[14]));
					
					    
					if(FieldNumbers[15]<0)
					    pt1 = NULL;
					else    
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[15]));
					    
					if(FieldNumbers[16]<0)
					    pt2 = NULL;
					else 
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[16]));

					
					returnptr = FieldNumbers;
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function BndMPoleSymplectic4E2Pass\n");
				}
		}


	irho = ba/le;
	
	BndMPoleSymplectic4E2Pass(r_in, le, irho, A, B, max_order, num_int_steps, 
								entrance_angle, exit_angle, fint1, fint2, gap, h1,h2,pt1, pt2, pr1, pr2, num_particles);
	

	
	return(returnptr);

}


 






void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	int m,n;
	double *r_in;
	double le, ba, *A, *B;  
	double irho;
	int max_order, num_int_steps;
	double entrance_angle, exit_angle ;
	double  *pr1, *pr2, *pt1, *pt2, fint1, fint2, gap,h1,h2;  
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
	/*else 		mexErrMsgTxt("Required field 'PolynomA' was not found in the element data structure");*/
				    
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

			    
		tmpmxptr = mxGetField(prhs[0],0,"H1");
	    if(tmpmxptr)
	        h1 = mxGetScalar(tmpmxptr);
	    else
	        h1 = 0;
		tmpmxptr = mxGetField(prhs[0],0,"H2");
	    
		if(tmpmxptr)
	        h2 = mxGetScalar(tmpmxptr);
	    else
	        h2 = 0;
	    
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
	BndMPoleSymplectic4E2Pass(r_in, le, irho, A, B, max_order, num_int_steps, 
								entrance_angle, exit_angle, fint1, fint2, gap, h1,h2,pt1, pt2, pr1, pr2, n);

	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(7,1);
	    
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("BendingAngle"));
	    mxSetCell(plhs[0],2,mxCreateString("EntranceAngle"));
	    mxSetCell(plhs[0],3,mxCreateString("ExitAngle"));
          /*mxSetCell(plhs[0],4,mxCreateString("PolynomA"));*/
	    mxSetCell(plhs[0],4,mxCreateString("PolynomB"));
	    mxSetCell(plhs[0],5,mxCreateString("MaxOrder"));
	    mxSetCell(plhs[0],6,mxCreateString("NumIntSteps"));	 	    
	    
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(9,1);
	        mxSetCell(plhs[1],0,mxCreateString("FullGap"));
	        mxSetCell(plhs[1],1,mxCreateString("FringeInt1"));
	        mxSetCell(plhs[1],2,mxCreateString("FringeInt2"));
	        mxSetCell(plhs[1],3,mxCreateString("H1"));
			mxSetCell(plhs[1],4,mxCreateString("H2"));
			mxSetCell(plhs[1],5,mxCreateString("T1"));
	        mxSetCell(plhs[1],6,mxCreateString("T2"));
	        mxSetCell(plhs[1],7,mxCreateString("R1"));
	        mxSetCell(plhs[1],8,mxCreateString("R2"));
						
	    }
	}



}
#endif /*MATLAB_MEX_FILE*/
