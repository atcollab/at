/* findmpoleraddifmatrix.c

   mex-function to calculate radiation diffusion matrix B defined in [2] 
   for multipole elements in MATLAB Accelerator Toolbox
   A.Terebilo 8/14/00

   References
   [1] M.Sands 'The Physics of Electron Storage Rings
   [2] Ohmi, Kirata, Oide 'From the beam-envelope matrix to synchrotron
   radiation integrals', Phys.Rev.E  Vol.49 p.751 (1994)
*/

#include "mex.h"
#include "atlalib.c"
#include <math.h>


/* Fourth order-symplectic integrator constants */

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

/* Physical constants used in the calculations */

#define TWOPI		6.28318530717959
#define CGAMMA 	8.846056192e-05 			/* [m]/[GeV^3] Ref[1] (4.1)      */
#define M0C2      5.10999060e5				/* Electron rest mass [eV]       */
#define LAMBDABAR 3.86159323e-13			/* Compton wavelength/2pi [m]    */
#define CER   		2.81794092e-15			/* Classical electron radius [m] */
#define CU        1.323094366892892			/* 55/(24*sqrt(3)) factor        */



#define SQR(X) ((X)*(X))



void smpledge(double* r, double inv_rho, double angle)
{	double psi = inv_rho*tan(angle);
	r[1]+=r[0]*psi;
	r[3]-=r[2]*psi;
}


double B2perp(double bx, double by, double irho, 
                            double x, double xpr, double y, double ypr)
/* Calculates sqr(|e x B|) , where e is a unit vector in the direction of velocity   */
    
{	double v_norm2;
	v_norm2 = 1/(SQR(1+x*irho)+ SQR(xpr) + SQR(ypr));

	/* components of the  velocity vector
	   double ex, ey, ez;
	   ex = xpr; 
	   ey = ypr; 
	   ez = (1+x*irho);
	*/
  	
	return((SQR(by*(1+x*irho)) + SQR(bx*(1+x*irho)) + SQR(bx*ypr - by*xpr) )*v_norm2) ;



} 
 

void thinkickrad(double* r, double* A, double* B, double L, double irho, double E0, int max_order)

/***************************************************************************** 
Calculate and apply a multipole kick to a phase space vector *r in a multipole element.
The reference coordinate system  may have the curvature given by the inverse 
(design) radius irho. irho = 0 for straight elements

IMPORTANT !!!
The desighn magnetic field Byo that provides this curvature By0 = irho * E0 /(c*e)
MUST NOT be included in the dipole term PolynomB(1)(MATLAB notation)(B[0] C notation) 
of the By field expansion
HOWEVER!!! to calculate the effect of classical radiation the full field must be 
used in the square of the |v x B|.
When calling B2perp(Bx, By, ...), use the By = ReSum + irho, where ReSum is the 
normalized vertical field - sum of the polynomial terms in PolynomB.

The kick is given by

           e L      L delta      L x
theta  = - --- B  + -------  -  -----  , 
     x     p    y     rho           2
            0                    rho

         e L
theta  = --- B
     y    p   x
           0

Note: in the US convention the field is written as:

                        max_order+1
                          ----
                          \                       n-1
	   (B + iB  ) = B rho  >   (ia  + b ) (x + iy)
        y    x           /       n    n
	                      ----
                          n=1

Use different index notation 

                         max_order
                           ----
                           \                       n
	   (B + iB  )/ B rho  =  >   (iA  + B ) (x + iy)
        y    x             /       n    n
  	                        ----
                           n=0

A,B: i=0 ... i=max_order
[0] - dipole, [1] - quadrupole, [2] - sextupole ...
units for A,B[i] = 1/[m]^(i+1)
Coeficients are stored in the PolynomA, PolynomB field of the element
structure in MATLAB


******************************************************************************/
{  int i;
 	double ImSum = A[max_order];
	double ReSum = B[max_order];	
	double x ,xpr, y, ypr, p_norm,dp_0, B2P;
	double ReSumTemp;
	double CRAD = CGAMMA*E0*E0*E0/(TWOPI*1e27);
  	
	/* recursively calculate the local transvrese magnetic field
	   Bx = ReSum, By = ImSum
	*/
	for(i=max_order-1;i>=0;i--)
		{	ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
			ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
			ReSum = ReSumTemp;
		}
	

	/* calculate angles from momentas */	
	p_norm = 1/(1+r[4]);
	x   = r[0];
	xpr = r[1]*p_norm;
	y   = r[2];
	ypr = r[3]*p_norm;


	B2P = B2perp(ImSum, ReSum +irho, irho, x , xpr, y ,ypr);
	
	dp_0 = r[4]; /* save a copy of the initial value of dp/p */

	r[4] = r[4] - CRAD*SQR(1+r[4])*B2P*(1 + x*irho + (SQR(xpr)+SQR(ypr))/2 )*L;

	/* recalculate momentums from angles after losing energy to radiation 	*/
	p_norm = 1/(1+r[4]);
	r[1] = xpr/p_norm;
	r[3] = ypr/p_norm;

  	
	r[1] -=  L*(ReSum-(dp_0-r[0]*irho)*irho);
	r[3] +=  L*ImSum;
	r[5] +=  L*irho*r[0]; /* pathlength */


}

void thinkickM(double* orbit_in, double* A, double* B, double L, 
							double irho, int max_order, double *M66)
/* Calculate the symplectic (no radiation) transfer matrix of a 
   thin multipole kick near the entrance point orbit_in
   For elements with straight coordinate system irho = 0
   For curved elements the B polynomial (PolynomB in MATLAB) 
   MUST NOT include the guide field  By0 = irho * E0 /(c*e)

   M is a (*double) pointer to a preallocated 1-dimentional array 
   of 36 elements of matrix M arranged column-by-column 
*/
{  int m,n;
   
	double ReSumNTemp;
   double ImSumN = max_order*A[max_order];
   double ReSumN = max_order*B[max_order];	

   /* Recursively calculate the derivatives
	   ReSumN = (irho/B0)*Re(d(By + iBx)/dx)
	   ImSumN = (irho/B0)*Im(d(By + iBx)/dy)
    */
	for(n=max_order-1;n>0;n--)
		{	ReSumNTemp = (ReSumN*orbit_in[0] - ImSumN*orbit_in[2]) + n*B[n];
			ImSumN = ImSumN*orbit_in[0] +  ReSumN*orbit_in[2] + n*A[n];
			ReSumN = ReSumNTemp;
		}

	/* Initialize M66 to a 6-by-6 identity matrix */
	for(m=0;m<36;m++)
		M66[m]= 0;
	/* Set diagonal elements to 1 */
	for(m=0;m<6;m++)
		M66[m*7] = 1;
	
	/* The relationship between indexes when a 6-by-6 matrix is 
	   represented in MATLAB as one-dimentional array containing
	   36 elements arranged column-by-column is
	   [i][j] <---> [i+6*j] 
	*/

	M66[1]   = -L*ReSumN;				/* [1][0] */
	M66[13]  =  L*ImSumN;				/* [1][2] */
	M66[3]   =  L*ImSumN;				/* [3][0] */
	M66[15]  =  L*ReSumN;				/* [3][2] */
	M66[25]  =  L*irho;					/* [1][4] */
	M66[1]  += -L*irho*irho;			/* [1][0] */
	M66[5]   =  L*irho;					/* [5][0] */

}



void thinkickB(double* orbit_in, double* A, double* B, double L, 
							double irho, int max_order, double E0, double *B66)

/* Calculate Ohmi's diffusion matrix of a thin multipole  element 
   For elements with straight coordinate system irho = 0
   For curved elements the B polynomial (PolynomB in MATLAB) 
   MUST NOT include the guide field  By0 = irho * E0 /(c*e)
   The result is stored in a preallocated 1-dimentional array B66
   of 36 elements of matrix B arranged column-by-column
*/

{	double BB,B2P,B3P;
	int i;
	double p_norm = 1/(1+orbit_in[4]);
	double p_norm2 = SQR(p_norm);
 	double ImSum = A[max_order];
	double ReSum = B[max_order];	
	double ReSumTemp;
	
  	/* recursively calculate the local transvrese magnetic field
	   ReSum = irho*By/B0
	   ImSum = irho*Bx/B0
	*/

	for(i=max_order-1;i>=0;i--)
		{	ReSumTemp = ReSum*orbit_in[0] - ImSum*orbit_in[2] + B[i];
			ImSum = ImSum*orbit_in[0] +  ReSum*orbit_in[2] + A[i];
			ReSum = ReSumTemp;
		}
	
	
	/* calculate |B x n|^3 - the third power of the B field component 
	   orthogonal to the normalized velocity vector n
    */
	B2P = B2perp(ImSum, ReSum +irho, irho, orbit_in[0] , orbit_in[1]*p_norm , 
								orbit_in[2] , orbit_in[3]*p_norm );
	B3P = B2P*sqrt(B2P);

	BB = CU * CER * LAMBDABAR *  pow(E0/M0C2,5) * L * B3P * SQR(SQR(1+orbit_in[4]))*
				(1+orbit_in[0]*irho + (SQR(orbit_in[1])+SQR(orbit_in[3]))*p_norm2/2);

	
	/* When a 6-by-6 matrix is represented in MATLAB as one-dimentional 
	   array containing 36 elements arranged column-by-column,
	   the relationship between indexes  is
	   [i][j] <---> [i+6*j] 

	*/
	
	/* initialize B66 to 0 */
	for(i=0;i<36;i++)
		B66[i] = 0;
	
	/* Populate B66 */
	B66[7]     = BB*SQR(orbit_in[1])*p_norm2;
	B66[19]    = BB*orbit_in[1]*orbit_in[3]*p_norm2;
	B66[9]     = BB*B66[19];
	B66[21]    = BB*SQR(orbit_in[3])*p_norm2;
	B66[10] = BB*orbit_in[1]*p_norm;
	B66[25] = BB*B66[10];
	B66[22] = BB*orbit_in[3]*p_norm;
	B66[27] = BB* B66[22];
	B66[28] = BB;
}





void drift_propagateB(double *orb_in, double L,  double *B)
{	/* Propagate cumulative Ohmi's diffusion matrix B through a drift
	   B is a (*double) pointer to 1-dimentional array 
	   containing 36 elements of matrix elements arranged column-by-column
	   as in MATLAB representation 

	   The relationship between indexes when a 6-by-6 matrix is 
	   represented in MATLAB as one-dimentional array containing
	   36 elements arranged column-by-column is
	   [i][j] <---> [i+6*j] 
	*/
		
	int m;
		
	double *DRIFTMAT = (double*)mxCalloc(36,sizeof(double));
	for(m=0;m<36;m++)
		DRIFTMAT[m] = 0;
	/* Set diagonal elements to 1	*/
	for(m=0;m<6;m++)
		DRIFTMAT[m*7] = 1;

	DRIFTMAT[6]  =  L/(1+orb_in[4]);
	DRIFTMAT[20] =  DRIFTMAT[6];
	DRIFTMAT[24] = -L*orb_in[1]/SQR(1+orb_in[4]);
	DRIFTMAT[26] = -L*orb_in[3]/SQR(1+orb_in[4]);
	DRIFTMAT[11] =  L*orb_in[1]/SQR(1+orb_in[4]);
	DRIFTMAT[23] =  L*orb_in[3]/SQR(1+orb_in[4]);	
	DRIFTMAT[29] = -L*(SQR(orb_in[1])+SQR(orb_in[3]))/((1+orb_in[4])*SQR(1+orb_in[4]));

	ATsandwichmmt(DRIFTMAT,B);
	mxFree(DRIFTMAT);
	
}


void edge_propagateB(double inv_rho, double angle, double *B)

{	/* Propagate  Ohmi's diffusion matrix B
	   through a focusing edge  B -> E*B*E'
	    where  E is a linear map of an edge 
	*/
	int m;
	double psi = inv_rho*tan(angle);
	
	for(m=0;m<6;m++)
		{	B[1+6*m] += psi*B[6*m];
			B[3+6*m] -= psi*B[2+6*m];
		}
	for(m=0;m<6;m++)
		{	B[m+6*1] += psi*B[m+6*0];
			B[m+6*3] -= psi*B[m+6*2];
		}
}

void FindElemB(double *orbit_in, double le, double irho, double *A, double *B,
					double *pt1, double* pt2,double *PR1, double *PR2,
					double entrance_angle, 	double exit_angle,	
					int max_order, int num_int_steps,
					double E0, double *BDIFF)

{	/* Find Ohmi's diffusion matrix BDIFF of a thick multipole
	   BDIFF - cumulative Ohmi's diffusion is initialized to 0
	   BDIFF is preallocated 1 dimensional array to store 6-by-6 matrix 
	   columnwise
	*/
	
	int m;	
	double  *MKICK, *BKICK;

	/* 4-th order symplectic integrator constants */
	double SL, L1, L2, K1, K2;
	SL = le/num_int_steps;
	L1 = SL*DRIFT1;
	L2 = SL*DRIFT2;
	K1 = SL*KICK1;
	K2 = SL*KICK2;
	
	
	/* Allocate memory for thin kick matrix MKICK
	   and a diffusion matrix BKICK
	*/
 	MKICK = (double*)mxCalloc(36,sizeof(double));
	BKICK = (double*)mxCalloc(36,sizeof(double));
	for(m=0; m < 6; m++)
		{	MKICK[m] = 0;
			BKICK[m] = 0;
		}
	
	/* Transform orbit to a local coordinate system of an element */
	
	ATaddvv(orbit_in,pt1);	
	ATmultmv(orbit_in,PR1);	

	/* This coordinate transformation does not affect 
	   the cumulative diffusion matrix BDIFF
	   E*BDIFF*E' :   BDIFF stays zero	

	*/	
	smpledge(orbit_in, irho, entrance_angle);	/* change in the input orbit 
												   from edge focusing
												*/
	
	edge_propagateB(irho,entrance_angle,BDIFF);		/* propagate the initial 
													   MRAD and BDIFF through 
													   the entrance edge
													*/

	/* Propagate orbit_in and BDIFF through a 4-th orderintegrator */

	for(m=0; m < num_int_steps; m++) /* Loop over slices	*/			
		{		drift_propagateB(orbit_in,L1, BDIFF);
				ATdrift6(orbit_in,L1);
				
				thinkickM(orbit_in, A,B, K1, irho, max_order, MKICK);
				thinkickB(orbit_in, A,B, K1, irho, max_order, E0, BKICK);
				ATsandwichmmt(MKICK,BDIFF);
				ATaddmm(BKICK,BDIFF);
				thinkickrad(orbit_in, A, B, K1, irho, E0, max_order);
		
				drift_propagateB(orbit_in,L2, BDIFF);
				ATdrift6(orbit_in,L2);
				
				thinkickM(orbit_in, A,B, K2, irho, max_order, MKICK);
				thinkickB(orbit_in, A,B, K2, irho, max_order, E0, BKICK);
				ATsandwichmmt(MKICK,BDIFF);
				ATaddmm(BKICK,BDIFF);
				thinkickrad(orbit_in, A, B, K2, irho, E0, max_order);
	
				drift_propagateB(orbit_in,L2, BDIFF);
				ATdrift6(orbit_in,L2);
				
				thinkickM(orbit_in, A,B, K1, irho, max_order, MKICK);
				thinkickB(orbit_in, A,B, K1, irho, max_order, E0, BKICK);
				ATsandwichmmt(MKICK,BDIFF);
				ATaddmm(BKICK,BDIFF);
				thinkickrad(orbit_in, A, B,  K1, irho, E0, max_order);

				drift_propagateB(orbit_in,L1, BDIFF);
				ATdrift6(orbit_in,L1);
		}  
		smpledge(orbit_in, irho, exit_angle);
		edge_propagateB(irho,exit_angle,BDIFF);
				
		ATsandwichmmt(PR2,BDIFF);
												
		mxFree(MKICK);
		mxFree(BKICK);
}


void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* The calling syntax of this mex-function from MATLAB is
   FindMPoleRadDiffMatrix(ELEMENT, ORBIT)
   ELEMENT is the element structure with field names consistent with 
           a multipole transverse field model.
   ORBIT is a 6-by-1 vector of the closed orbit at the entrance (calculated elsewhere)
*/
{	int m,n;
	double le, ba, *A, *B;  
	double irho;
	const mxArray * globvalptr;
	mxArray *E0ptr;
	double E0;		/* Design energy [eV] to be obtained from MATLAB global workspace */
	int max_order, num_int_steps;
	double entrance_angle, exit_angle ;
	double *BDIFF;
	mxArray  *mxtemp;

	double *orb, *orb0;
	double *pt1, *pt2, *PR1, *PR2;


	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(!(m==6 && n==1))
		mexErrMsgTxt("Second argument must be a 6-by-1 column vector");
    
	/* ALLOCATE memory for the output array */
	plhs[0] = mxCreateDoubleMatrix(6,6,mxREAL);
	BDIFF = mxGetPr(plhs[0]);


	/* If the ELEMENT sructure does not have fields PolynomA and PolynomB
	   return zero matrix and  exit
	*/
	if(mxGetField(prhs[0],0,"PolynomA") == NULL ||  mxGetField(prhs[0],0,"PolynomB") == NULL)
		return;
	

	/* retrieve the value of design Energy [GeV]
	   contained in MATLAB global variable GLOBVAL.
	   GLOBVAL is a MATLAB structure
	   GLOBVAL.E0 contains the design energy of the ring [eV]
    */

	globvalptr=mexGetArrayPtr("GLOBVAL","global");
	if(globvalptr != NULL)
		{	E0ptr = mxGetField(globvalptr,0,"E0");
			if(E0ptr !=NULL)
				E0 = mxGetScalar(E0ptr);
			else
				mexErrMsgTxt("Global variable GLOBVAL does not have a field 'E0'");
		}
	else
		mexErrMsgTxt("Global variable GLOBVAL does not exist");

	orb0 = mxGetPr(prhs[1]);
	/* make local copy of the input closed orbit vector */
	orb = (double*)mxCalloc(6,sizeof(double));
	for(m=0;m<6;m++)
		orb[m] = orb0[m];
    
	/* Retrieve element information */
	
	le = mxGetScalar(mxGetField(prhs[0],0,"Length"));
	
	/* If ELEMENT has a zero length, return zeros matrix end exit */
	if(le == 0)
		return;
	
	A = mxGetPr(mxGetField(prhs[0],0,"PolynomA"));
	B = mxGetPr(mxGetField(prhs[0],0,"PolynomB"));

	

		
	mxtemp = mxGetField(prhs[0],0,"NumIntSteps");
   if(mxtemp != NULL)
		num_int_steps = (int)mxGetScalar(mxtemp);
	else
		mexErrMsgTxt("Field 'NumIntSteps' not found in the ELEMENT structure");

	mxtemp = mxGetField(prhs[0],0,"MaxOrder");
   if(mxtemp != NULL)
		max_order = (int)mxGetScalar(mxtemp);
	else
		mexErrMsgTxt("Field 'MaxOrder' not found in the ELEMENT structure");


	mxtemp = mxGetField(prhs[0],0,"BendingAngle");
   if(mxtemp != NULL)
		{	ba = mxGetScalar(mxtemp);
			irho = ba/le;
		}
	else
		{	ba = 0;
			irho = 0;
		}
		
	mxtemp = mxGetField(prhs[0],0,"EntranceAngle");
	if(mxtemp != NULL)
		entrance_angle = mxGetScalar(mxtemp);
	else
			entrance_angle =0;

	mxtemp = mxGetField(prhs[0],0,"ExitAngle");
	if(mxtemp != NULL)
		exit_angle = mxGetScalar(mxtemp);
	else
			exit_angle =0;

	pt1 = mxGetPr(mxGetField(prhs[0],0,"T1"));
	pt2 = mxGetPr(mxGetField(prhs[0],0,"T2"));
	PR1 = mxGetPr(mxGetField(prhs[0],0,"R1"));
	PR2 = mxGetPr(mxGetField(prhs[0],0,"R2"));
	

	FindElemB(orb, le, irho, A, B, 
					pt1, pt2, PR1, PR2,
					entrance_angle, 	exit_angle, 
					max_order, num_int_steps, E0, BDIFF);
}


