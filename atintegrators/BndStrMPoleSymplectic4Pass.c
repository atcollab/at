#include "atelem.c"
#include "atlalib.c"
#include "atphyslib.c"

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

/* Straight dipole w/ multipole using Symplectic Integration and rotation at
 * dipole faces.
 * Created by Xiaobiao Huang, 7/31/2018 */
#define SQR(X) ((X)*(X))

/* At entrance edge: move particles to the field edge and convert coordinates
 * to x, dx/dz, y, dy/dz, then convert to x, px, y, py as integration is done
 * with px, py */
void E1rotation(double *r,double X0ref, double E1)
{
    double x0,dxdz0, dydz0, psi;
    double fac;
    
    dxdz0 = r[1]/sqrt(SQR(1+r[4])-SQR(r[1])-SQR(r[3]));
    dydz0 = r[3]/sqrt(SQR(1+r[4])-SQR(r[1])-SQR(r[3]));
    x0 = r[0];
    
    psi = atan(dxdz0);
    r[0] = r[0]*cos(psi)/cos(E1+psi)+X0ref;
    r[1] = tan(E1+psi);
    r[3] = dydz0/(cos(E1)-dxdz0*sin(E1));
    r[2] += x0*sin(E1)*r[3];
    
    r[5] += x0*tan(E1)/(1-dxdz0*tan(E1))*sqrt(1+SQR(dxdz0)+SQR(dydz0));
    
    /* convert to px, py */
    fac = sqrt(1+SQR(r[1])+SQR(r[3]));
    r[1] = r[1]*(1+r[4])/fac;
    r[3] = r[3]*(1+r[4])/fac;
}
/* At Exit Edge: : move particles to arc edge and convert coordinates to x, px,
 * y, py */
void E2rotation(double *r,double X0ref, double E2)
{
    double x0;
    double dxdz0, dydz0, psi, fac;
    
    dxdz0 = r[1]/sqrt(SQR(1+r[4])-SQR(r[1])-SQR(r[3]));
    dydz0 = r[3]/sqrt(SQR(1+r[4])-SQR(r[1])-SQR(r[3]));
    x0 = r[0];
    
    psi = atan(dxdz0);
    fac = sqrt(1+SQR(dxdz0)+SQR(dydz0));
    
    r[0] = (r[0]-X0ref)*cos(psi)/cos(E2+psi);
    r[1] = tan(E2+psi);
    r[3] = dydz0/(cos(E2)-dxdz0*sin(E2));
    r[2] += r[3]*(x0-X0ref)*sin(E2);
    
    r[5] += (x0-X0ref)*tan(E2)/(1-dxdz0*tan(E2))*fac;

    /* convert to px, py */
    fac = sqrt(1+SQR(r[1])+SQR(r[3]));
    r[1] = r[1]*(1+r[4])/fac;
    r[3] = r[3]*(1+r[4])/fac;
    
}

void edgey(double* r, double inv_rho, double edge_angle)
{	/* Edge focusing in dipoles with hard-edge field for vertical only */
    double psi = inv_rho*tan(edge_angle);

    /* r[1]+=r[0]*psi; */
	r[3]-=r[2]*psi;
}


void edgey_fringe(double* r, double inv_rho, double edge_angle, double fint, double gap)
{   /* Edge focusing in dipoles with fringe field, for vertical only */
    double fx = inv_rho*tan(edge_angle);
    double psi_bar = edge_angle-inv_rho*gap*fint*(1+sin(edge_angle)*sin(edge_angle))/cos(edge_angle)/(1+r[4]);
    double fy = inv_rho*tan(psi_bar);
    /* r[1]+=r[0]*fx; */
    r[3]-=r[2]*fy;    
}


void ladrift6(double* r, double L)
/* large angle drift, X. Huang, 7/31/2018  
 * Input parameter L is the physical length
 * 1/(1+delta) normalization is done internally
 * Hamiltonian H = (1+\delta)-sqrt{(1+\delta)^2-p_x^2-p_y^2}, change sign for
 * $\Delta z$ in AT */
{	double p_norm = 1./sqrt(SQR(1+r[4])-SQR(r[1])-SQR(r[3])); 
	double NormL  = L*p_norm;   
	r[0]+= NormL*r[1]; 
	r[2]+= NormL*r[3];
	r[5]+= L*(p_norm*(1+r[4])-1.);
}


void bndstrthinkick(double* r, double* A, double* B, double L, double irho, int max_order)

/***************************************************************************** 
Calculate multipole kick in a straight bending magnet, This is not the usual Bends!
 *created by X. Huang, 7/31/2018
The reference coordinate system  is straight in s. 
 *The B vector does not contain b0, we assume b0=irho

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
	   Bx = ReSum, By = ImSum */
    B[0] = irho; 
	for(i=max_order-1;i>=0;i--)
		{	ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
			ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
			ReSum = ReSumTemp;
		}
	
  	
	r[1] -=  L*(ReSum);
	r[3] +=  L*ImSum;
	r[5] +=  0; /* pathlength */


}



void BndStrMPoleSymplectic4Pass(double *r, double le, double irho, double *A, double *B,
					int max_order, int num_int_steps,
					double entrance_angle, 	double exit_angle,double X0ref, double ByError,double RefDZ,
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
		
    /* mexPrintf("E0ref=%f\n",X0ref); */
    
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
			            edgey_fringe(r6, irho+B[1]*X0ref, entrance_angle,fint1,gap);
			        else
			            edgey(r6, irho+B[1]*X0ref, entrance_angle);
				 	
                    /* Rotate and translate to straight Cartesian coordinate */
                    E1rotation(r6, X0ref, entrance_angle);
                    
                    /* integrator */
					for(m=0; m < num_int_steps; m++) /* Loop over slices */
						{		r6 = r+c*6;	
								
								ladrift6(r6,L1);
           					    bndstrthinkick(r6, A, B, K1, irho, max_order);
								ladrift6(r6,L2);
           					    bndstrthinkick(r6, A, B, K2, irho, max_order);
								ladrift6(r6,L2);
		     					bndstrthinkick(r6, A, B,  K1, irho, max_order);
								ladrift6(r6,L1);	
						}  
					
                     /* Rotate and translate back to curvilinear coordinate */
                    E2rotation(r6, X0ref, exit_angle);
                    r6[5] -= RefDZ;
                    
					if(useFringe2)
			            edgey_fringe(r6, irho+B[1]*X0ref, exit_angle,fint2,gap);
			        else
			            edgey(r6, irho+B[1]*X0ref, exit_angle);	
					/* edge focus */


                   
                    
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
    double X0ref, ByError,flen,RefDZ;

	int max_order, num_int_steps;
	double le,ba,irho;
	int *returnptr;
	int *NewFieldNumbers, fnum;

	ByError =0;
    
	switch(mode)
		{   case MAKE_LOCAL_COPY: 	/* Find field numbers first, save a list
                                       of field number in an array and make
                                       returnptr point to that array */
				{	
					/* Allocate memory for integer array of field numbers for
                       faster future reference */
		
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
					
				
                    fnum = mxGetFieldNumber(ElemData,"X0ref");
					NewFieldNumbers[15] = fnum;
					if(fnum<0) 
					    X0ref = 0;
					else
					    X0ref = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
                    
                    fnum = mxGetFieldNumber(ElemData,"ByError");
					NewFieldNumbers[16] = fnum;
					if(fnum<0) 
					    ByError = 0;
					else
					    ByError = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
                    
                    fnum = mxGetFieldNumber(ElemData,"RefDZ");
					NewFieldNumbers[17] = fnum;
					if(fnum<0) 
					    RefDZ = 0;
					else
					    RefDZ = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
                    
					returnptr = NewFieldNumbers;

				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field
                                       numbers, the second argument ponter to
                                       the array of field numbers is previously
                                       created with:
                                       QuadLinPass( ..., MAKE_LOCAL_COPY) */
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
					
                    if(FieldNumbers[15]<0) 
					    X0ref = 0;
					else
					    X0ref = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[15]));
                    
					if(FieldNumbers[16]<0) 
					    ByError = 0;
					else
					    ByError = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[16]));
                    
                    if(FieldNumbers[17]<0) 
					    RefDZ = 0;
					else
					    RefDZ = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[17]));
                    
					returnptr = FieldNumbers;
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function BndStrMPoleSymplectic4Pass\n");
				}
		}


	irho = ba/le;
	flen = 2.0/irho*sin(ba/2.0); /* field length */
	BndStrMPoleSymplectic4Pass(r_in, flen, irho, A, B, max_order, num_int_steps, 
								entrance_angle, exit_angle, X0ref, ByError,RefDZ,fint1, fint2, gap, pt1, pt2, pr1, pr2, num_particles);
	

	
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
    double X0ref, ByError, RefDZ;
    double flen;
    
    mxArray *tmpmxptr;

    X0ref = 0;
    ByError = 0;
    
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
		
		tmpmxptr = mxGetField(prhs[0],0,"X0ref");
	    if(tmpmxptr)
	        X0ref = mxGetScalar(tmpmxptr);
	    else
	        X0ref = 0;
        
        tmpmxptr = mxGetField(prhs[0],0,"ByError");
	    if(tmpmxptr)
	        ByError = mxGetScalar(tmpmxptr);
	    else
	        ByError = 0;
 
        tmpmxptr = mxGetField(prhs[0],0,"RefDZ");
	    if(tmpmxptr)
	        RefDZ = mxGetScalar(tmpmxptr);
	    else
	        RefDZ = 0;
        
    irho = ba/le;
    flen = 2.0/irho*sin(ba/2.0); /* field length */
    
    plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	BndStrMPoleSymplectic4Pass(r_in, flen, irho, A, B, max_order, num_int_steps, 
								entrance_angle, exit_angle, X0ref, ByError, RefDZ, fint1, fint2, gap, pt1, pt2, pr1, pr2, n);

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
            mxSetCell(plhs[1],7,mxCreateString("X0ref"));
	        mxSetCell(plhs[1],8,mxCreateString("ByError"));
            mxSetCell(plhs[1],9,mxCreateString("RefDZ"));
	    }
	}



}



