#include "atelem.c"
#include "atlalib.c"
#include "atphyslib.c"

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656


struct elem
{
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
    double h1;
    double h2;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
};


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
{
    int i;
    double ReSum = B[max_order];
    double ImSum = A[max_order];
    double ReSumTemp;
  	/* recursively calculate the local transvrese magnetic field
	   Bx = ReSum, By = ImSum
	*/
    for(i=max_order-1;i>=0;i--)
    {
        ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
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
        double *R1, double *R2, 
        double *RApertures, double *EApertures, int num_particles)
{
    int c,m;
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
    
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(shared) shared(r,num_particles) private(c,r6,m)
    for(c = 0;c<num_particles;c++)	/* Loop over particles  */
    {
        r6 = r+c*6;
        if(!atIsNaN(r6[0]))
        {
            /*  misalignment at entrance  */
            if(useT1) ATaddvv(r6,T1);
            if(useR1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* edge focus */
            if(useFringe1)
            {
                edge_fringe2A(r6, irho, entrance_angle,fint1,gap,h1,B[1]);
            }
            else
            {
                edge_fringe2A(r6, irho, entrance_angle,0,0,h1,B[1]);
            }
            /* integrator */
            for(m=0; m < num_int_steps; m++) /* Loop over slices*/
            {
                r6 = r+c*6;
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
            }
            else
            {
                edge_fringe2B(r6, irho, exit_angle,0,0,h2,B[1]);
            }
            /* Check physical apertures at the exit of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* Misalignment at exit */
            if(useR2) ATmultmv(r6,R2);
            if(useT2) ATaddvv(r6,T2);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double irho;
    if (!Elem) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, Gap, Fint1, Fint2;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, h1, h2, *R1, *R2, *T1, *T2, *EApertures, *RApertures;
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
        h1=atGetOptionalDouble(ElemData,"H1", 0); check_error();
        h2=atGetOptionalDouble(ElemData,"H2", 0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
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
        Elem->h1=h1;
        Elem->h2=h2;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
    }
    irho = Elem->BendingAngle / Elem->Length;
    BndMPoleSymplectic4E2Pass(r_in, Elem->Length, irho, Elem->PolynomA, 
            Elem->PolynomB, Elem->MaxOrder, Elem->NumIntSteps,
            Elem->EntranceAngle, Elem->ExitAngle, Elem->FringeInt1,
            Elem->FringeInt2, Elem->FullGap, Elem->h1, Elem->h2,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2, 
            Elem->RApertures,Elem->EApertures,num_particles);
    return Elem;
}

MODULE_DEF(BndMPoleSymplectic4E2Pass)        /* Dummy module initialisation */

#endif /*MATLAB_MEX_FILE || PYAT*/

#ifdef MATLAB_MEX_FILE

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double irho;
        double Length, BendingAngle, EntranceAngle, ExitAngle, Gap, Fint1, Fint2;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, h1, h2, *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgTxt("Second argument must be a 6 x N matrix");
        
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
        h1=atGetOptionalDouble(ElemData,"H1", 0); check_error();
        h2=atGetOptionalDouble(ElemData,"H2", 0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        irho = BendingAngle/Length;
        
        /* ALLOCATE memory for the output array of the same size as the input */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        BndMPoleSymplectic4E2Pass(r_in, Length, irho, PolynomA, PolynomB, MaxOrder, NumIntSteps, EntranceAngle, ExitAngle,
                Fint1, Fint2, Gap, h1, h2, T1, T2, R1, R2, RApertures,EApertures,num_particles);
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
        if (nlhs > 1) {    /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(11,1);
            mxSetCell(plhs[1],0,mxCreateString("FullGap"));
            mxSetCell(plhs[1],1,mxCreateString("FringeInt1"));
            mxSetCell(plhs[1],2,mxCreateString("FringeInt2"));
            mxSetCell(plhs[1],3,mxCreateString("H1"));
            mxSetCell(plhs[1],4,mxCreateString("H2"));
            mxSetCell(plhs[1],5,mxCreateString("T1"));
            mxSetCell(plhs[1],6,mxCreateString("T2"));
            mxSetCell(plhs[1],7,mxCreateString("R1"));
            mxSetCell(plhs[1],8,mxCreateString("R2"));
            mxSetCell(plhs[1],9,mxCreateString("RApertures"));
            mxSetCell(plhs[1],10,mxCreateString("EApertures"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/
