function [B, M, O] = findelemraddifmat(ELEM,orbit,varargin)
%FINDELEMRADDIFMAT calculates element 'radiation diffusion matrix' B
% [B, M, ORBITOUT] = FINDELEMRADDIFMAT(ELEM, ORBITIN);
% Ohmi, Kirata, Oide 'From the beam-envelope matrix to synchrotron
%    radiation integrals', Phys.Rev.E  Vol.49 p.751 (1994)


% Fourth order-symplectic integrator constants

DRIFT1   = 0.6756035959798286638
DRIFT2   = -0.1756035959798286639
KICK1    =  1.351207191959657328
KICK2    = -1.702414383919314656

% Physical constants used in calculations

TWOPI       = 6.28318530717959
CGAMMA      = 8.846056192e-05 			% [m]/[GeV^3] Ref[1] (4.1)  
M0C2        = 5.10999060e5              % Electron rest mass [eV]
LAMBDABAR   = 3.86159323e-13			% Compton wavelength/2pi [m]
CER         = 2.81794092e-15            % Classical electron radius [m]
CU          = 1.323094366892892         % 55/(24*sqrt(3))


function b2 = B2perp(B, irho, r6)
% Calculates sqr(|e x B|) , where e is a unit vector in the direction of
% velocity. Components of the  velocity vector:
% ex = xpr; 
% ey = ypr; 
% ez = (1+x*irho);

{	E = [r(2)/(1+r(5));r(4)/(1+r(5));1+r(1)*irho];
    b2 = sum(cross(E/norm(E),B).^2);
    b2 = dot(c,c);
} 
 

function rout = thinkickrad(rin, PolynomA, PolynomB, L, irho, E0, max_order)
% Propagate particle through a thin multipole with radiation
% Calculate field from polynomial coefficients
P = i*PolynomA(1:max_order+1)+PolynomB(1:max_order+1);
Z = cumprod([1, (rin(1)+i*rin(3))*ones(1,max_order)]);
S = sum(P.*Z);
Bx = real(S); By = imag(S);

B2P = B2perp([Bx By +irho 0], irho, r);
CRAD = CGAMMA*ELEM.Energy^3/(TWOPI*1e27);

% Propagate particle
rout = rin;

% Loss of energy (dp/p) due to radiation
rout(5) = rin(5) - CRAD*(1+rin(5))^2*B2P*...
    (1+rin(1)*irho + (rin(1)^2+rin(3)^2)/2/(1+rin(5))^2)*L;

% Change in transverse momentum due to radiation
%   Angle does not change but dp/p changes due to radiation
%   and therefore transverse canonical momentum changes 
%   px = x'*(1+dp/p)
%   py = y'*(1+dp/p)
rout(2 4]) = rin([2 4])*(1+rout(5))/(1+rin(5));

% transverse kick due to magnetic field
rout(2) = rout(2) - L*(Bx-(rin(5)-rin(1)*irho)*irho);
rout(4) = rout(4) + L*By;

% pathlength
rout(6) = rout(6) + L*irho*rin(1); 



function M = thinkickM(rin, PolynomA, PolynomB, L, irho, max_order)
%     Calculate the symplectic (no radiation) transfer matrix of a 
%    thin multipole kick near the entrance point orbit_in
%    For elements with straight coordinate system irho = 0
%    For curved elements the B polynomial (PolynomB in MATLAB) 
%    MUST NOT include the guide field  By0 = irho * E0 /(c*e)

{   P = i*PolynomA(2:max_order+1)+PolynomB(2:max_order+1);
    Z = cumprod([1, (rin(1)+i*rin(3))*ones(1,max_order-1)]);
    dB = sum(P.*(1:max_order).*Z);

    M = eye(6);



	M(2,1)   = -L*real(dB);
	M(2,3)   =  L*imag(dB);
	M(4,1)   =  L*imag(dB);
	M(4,3)   =  L*real(dB);
	M(2,5)   =  L*irho;
	M(2,1)   =  M(2,1) - L*irho*irho;
	M(6,1)   =  L*irho;

}



function B66 = thinkickB(rin, PolynomA, PolynomB, L, irho, E0, max_order)
%    Calculate Ohmi's diffusion matrix of a thin multipole  element 
%    For elements with straight coordinate system irho = 0
%    For curved elements the B polynomial (PolynomB in MATLAB) 
%    MUST NOT include the guide field  By0 = irho * E0 /(c*e)
%    The result is stored in a preallocated 1-dimentional array B66
%    of 36 elements of matrix B arranged column-by-column

P = i*PolynomA(1:max_order+1)+PolynomB(1:max_order+1);
Z = cumprod([1, (rin(1)+i*rin(3))*ones(1,max_order)]);
S = sum(P.*Z);
Bx = real(S); By = imag(S);

B2P = B2perp([Bx By +irho 0], irho, r);
B3P = B2P^(3/2);

p_norm = 1/(1+rin(5));
p_norm2 = p_norm^2;

BB = CU * CER * LAMBDABAR *  pow(E0/M0C2,5) * L * B3P * (1+rin(5))^4*
				(1+rin(1)*irho + (rin(2)^2+rin(4)^2)*p_norm2/2);


B66 = zeros(6);
B66(2,2)    = BB*rin(2)^2*p_norm2;
B66(2,4)    = BB*rin(2)*rin(4)*p_norm2;
B66(4,2)    = B66(2,4);
B66(4,4)    = BB*rin(4)^2*p_norm2;
B66(5,2)    = BB*rin(2)*p_norm;
B66(2,5)    = B66(5,2);
B66(5,4)    = BB*rin(4)*p_norm;
B66(4,5)    = B66(5,4);
B66(5,5)    = BB;


function = mvoid drift_propagateB(double *orb_in, double L,  double *B)
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


