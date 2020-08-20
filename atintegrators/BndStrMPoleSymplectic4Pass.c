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

struct elem {
    double Length;
    double BendingAngle;
    double EntranceAngle;
    double ExitAngle;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    /* Optional fields */
    double FullGap;
    double FringeInt1;
    double FringeInt2;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double X0ref;
    double ByError;
    double RefDZ;
};

void E1rotation(double *r,double X0ref, double E1)
/* At Entrance Edge:
 * move particles to the field edge and convert coordinates to x, dx/dz, y,
 * dy/dz, then convert to x, px, y, py as integration is done with px, py */
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

void E2rotation(double *r,double X0ref, double E2)
/* At Exit Edge:
 * move particles to arc edge and convert coordinates to x, px, y, py */
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
/* Edge focusing in dipoles with hard-edge field for vertical only */
{
    double psi = inv_rho*tan(edge_angle);
    /*r[1]+=r[0]*psi;*/
    r[3]-=r[2]*psi;
}

void edgey_fringe(double* r, double inv_rho, double edge_angle, double fint, double gap)
/* Edge focusing in dipoles with fringe field, for vertical only */
{
    double fx = inv_rho*tan(edge_angle);
    double psi_bar = edge_angle-inv_rho*gap*fint*(1+sin(edge_angle)*sin(edge_angle))/cos(edge_angle)/(1+r[4]);
    double fy = inv_rho*tan(psi_bar);
    /*r[1]+=r[0]*fx;*/
    r[3]-=r[2]*fy;
}

void ladrift6(double* r, double L)
/* large angle drift, X. Huang, 7/31/2018
 * Input parameter L is the physical length
 * 1/(1+delta) normalization is done internally
 * Hamiltonian H = (1+\delta)-sqrt{(1+\delta)^2-p_x^2-p_y^2}, change sign for
 * $\Delta z$ in AT */
{
    double p_norm = 1./sqrt(SQR(1+r[4])-SQR(r[1])-SQR(r[3]));
    double NormL = L*p_norm;
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
{
    int i;
	double ReSum = B[max_order];
 	double ImSum = A[max_order];
	double ReSumTemp;

  	/* recursively calculate the local transvrese magnetic field
	   Bx = ReSum, By = ImSum */
    B[0] = irho;
	for(i=max_order-1;i>=0;i--) {
	    ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
	    ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
	    ReSum = ReSumTemp;
    }
	r[1] -=  L*(ReSum);
	r[3] +=  L*ImSum;
	r[5] +=  0; /* pathlength */
}

void BndStrMPoleSymplectic4Pass(double *r, double le, double irho, double *A,
            double *B, int max_order, int num_int_steps,
            double entrance_angle, double exit_angle, double X0ref,
            double ByError, double RefDZ, double fint1, double fint2,
            double gap, double *T1, double *T2, double *R1, double *R2,
            int num_particles)
{
	int c, m;
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

	for(c = 0;c<num_particles;c++)	{   /* Loop over particles  */
        r6 = r+c*6;
	    if(!atIsNaN(r6[0])) {
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
			for(m=0; m < num_int_steps; m++) {  /* Loop over slices */
                r6 = r+c*6;
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
	        else    /* edge focus */
	            edgey(r6, irho+B[1]*X0ref, exit_angle);
			/* Misalignment at exit */
	        if(useR2)
	            ATmultmv(r6,R2);
            if(useT2)
	            ATaddvv(r6,T2);
        }
    }
}

/********** END PHYSICS SECTION **********************************************/
/*****************************************************************************/

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double irho, flen;
    if (!Elem) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, Gap, Fint1, Fint2, X0ref, ByError, RefDZ;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2;
        Length=atGetDouble(ElemData,"Length"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        /*optional fields*/
        Gap=atGetOptionalDouble(ElemData,"FullGap", 0); check_error();
        Fint1=atGetOptionalDouble(ElemData,"FringeInt1", 0); check_error();
        Fint2=atGetOptionalDouble(ElemData,"FringeInt2", 0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        X0ref=atGetOptionalDouble(ElemData,"X0ref", 0); check_error();
        ByError=atGetOptionalDouble(ElemData,"ByError", 0); check_error();
        RefDZ=atGetOptionalDouble(ElemData,"RefDZ", 0); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->BendingAngle=BendingAngle;
        Elem->EntranceAngle=EntranceAngle;
        Elem->ExitAngle=ExitAngle;
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        Elem->NumIntSteps=NumIntSteps;
        /*optional fields*/
        Elem->FullGap=Gap;
        Elem->FringeInt1=Fint1;
        Elem->FringeInt2=Fint2;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->X0ref=X0ref;
        Elem->ByError=ByError;
        Elem->RefDZ=RefDZ;
    }
    irho = Elem->BendingAngle / Elem->Length;
    flen = 2.0 / irho * sin(Elem->BendingAngle/2.0);
    BndStrMPoleSymplectic4Pass(r_in, flen, irho, Elem->PolynomA,
            Elem->PolynomB, Elem->MaxOrder, Elem->NumIntSteps,
            Elem->EntranceAngle, Elem->ExitAngle, Elem->X0ref, Elem->ByError,
            Elem->RefDZ, Elem->FringeInt1, Elem->FringeInt2, Elem->FullGap,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2, num_particles);
    return Elem;
}

MODULE_DEF(BndStrMPoleSymplectic4Pass)      /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double irho, flen;
        double Length, BendingAngle, EntranceAngle, ExitAngle, Gap, Fint1, Fint2, X0ref, ByError, RefDZ;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        Length=atGetDouble(ElemData,"Length"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        /*optional fields*/
        Gap=atGetOptionalDouble(ElemData,"FullGap", 0); check_error();
        Fint1=atGetOptionalDouble(ElemData,"FringeInt1", 0); check_error();
        Fint2=atGetOptionalDouble(ElemData,"FringeInt2", 0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        X0ref=atGetOptionalDouble(ElemData,"X0ref", 0); check_error();
        ByError=atGetOptionalDouble(ElemData,"ByError", 0); check_error();
        RefDZ=atGetOptionalDouble(ElemData,"RefDZ", 0); check_error();
        irho = BendingAngle/Length;
        flen = 2.0/irho*sin(BendingAngle/2.0); /* field length */
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        BndStrMPoleSymplectic4Pass(r_in, flen, irho, PolynomA, PolynomB,
            MaxOrder, NumIntSteps, EntranceAngle, ExitAngle, X0ref, ByError,
            RefDZ, Fint1, Fint2, Gap, T1, T2, R1, R2, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
	    plhs[0] = mxCreateCellMatrix(8,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("BendingAngle"));
	    mxSetCell(plhs[0],2,mxCreateString("EntranceAngle"));
	    mxSetCell(plhs[0],3,mxCreateString("ExitAngle"));
        mxSetCell(plhs[0],4,mxCreateString("PolynomA"));
	    mxSetCell(plhs[0],5,mxCreateString("PolynomB"));
	    mxSetCell(plhs[0],6,mxCreateString("MaxOrder"));
	    mxSetCell(plhs[0],7,mxCreateString("NumIntSteps"));
	    if(nlhs>1) {    /* list of optional fields */
	        plhs[1] = mxCreateCellMatrix(7,1);
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
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/
