#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "atphyslib.c"

#define SQR(X) ((X)*(X))

struct elem
{
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    double BendingAngle;
    double EntranceAngle;
    double ExitAngle;
    double Energy;
    /* Optional fields */
    double FringeInt1;
    double FringeInt2;
    double FullGap;
    double Scaling;
    double h1;
    double h2;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *KickAngle;
};

/*
 This code was modified from the original BndMPoleSymplectic4RadPass.c of AT to correctly integrate the Hamiltonian in 
 the curvilinear coordinate system of the dipole and to include the second order Transport map of the fringe field. Also 
 modified is the field Bx, By to include the curvature effect. 
 New version created by Xiaobiao Huang on 08/13/2009.
 Last modified on 8/26/2009
 */

static double B2perp(double bx, double by, double irho, 
                            double x, double xpr, double y, double ypr)
/* Calculates sqr(|e x B|) , where e is a unit vector in the direction of velocity  */    
{
    double v_norm2;
    v_norm2 = 1/(SQR(1+x*irho)+ SQR(xpr) + SQR(ypr));

	/* components of the  velocity vector
	   double ex, ey, ez;
	    ex = xpr; 
	    ey = ypr; 
	    ez = (1+x*irho);
	*/
    return((SQR(by*(1+x*irho)) + SQR(bx*(1+x*irho)) + SQR(bx*ypr - by*xpr) )*v_norm2) ;
}
 

static void bndthinkickrad(double* r, double* A, double* B, double L, double h, double E0,int max_order)
/*****************************************************************************
(1) PolynomA is neglected.
(2) The vector potential is expanded up to 4th order of x and y. 
(3) Coefficients in PolynomB higher than 4th order is treated as if they are on straight geometry.
*/
{
    int i;
    double ReSum = 0; /*B[max_order];*/
    double ImSum = 0; /*A[max_order];*/
    
    double ReSumTemp;
    double K1,K2;
    double x ,xpr, y, ypr, p_norm, B2P;
    
    double CRAD = CGAMMA*E0*E0*E0/(TWOPI*1e27);	/* [m]/[GeV^3] M.Sands (4.1) */
    
    K1 = B[1];
    K2 = (max_order>=2) ? B[2] : 0;
    
    ReSum = B[max_order];
    for(i=max_order-1;i>=0;i--) {
        ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
        ImSum = ImSum*r[0] +  ReSum*r[2] ;
        ReSum = ReSumTemp;
    }
    /* calculate angles from momentums 	*/
    p_norm = 1/(1+r[4]);
    x   = r[0];
    xpr = r[1]*p_norm;
    y   = r[2];
    ypr = r[3]*p_norm;
    /* see Iselin Part. Accel. 1985  */
    ImSum += h*(K1*h-K2)*y*y*y/6.0;
    ReSum += -K1*h*y*y/2.0 + h*(K1*h-K2)*x*y*y/2.0;
    
    B2P = B2perp(ImSum, ReSum +h, h, x , xpr, y ,ypr);
    
    r[4] = r[4] - CRAD*SQR(1+r[4])*B2P*(1 + x*h + (SQR(xpr)+SQR(ypr))/2 )*L;
    
    /* recalculate momentums from angles after losing energy for radiation 	*/
    p_norm = 1/(1+r[4]);
    r[1] = xpr/p_norm;
    r[3] = ypr/p_norm;
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

void BndMPoleSymplectic4E2RadPass(double *r, double le, double irho, double *A, double *B,
        int max_order, int num_int_steps,
        double entrance_angle, 	double exit_angle,
        double fint1, double fint2, double gap,double h1,double h2,
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double *KickAngle, double scaling, double E0, int num_particles)
{
    double SL = le/num_int_steps;
    double L1 = SL*DRIFT1;
    double L2 = SL*DRIFT2;
    double K1 = SL*KICK1;
    double K2 = SL*KICK2;
    bool useFringe1 = (fint1 != 0) && (gap != 0);
    bool useFringe2 = (fint2 != 0) && (gap != 0);
    double B0 = B[0];
    double A0 = A[0];

    if (KickAngle) {   /* Convert corrector component to polynomial coefficients */
        B[0] -= sin(KickAngle[0])/le;
        A[0] += sin(KickAngle[1])/le;
    }
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
    shared(r,num_particles,R1,T1,R2,T2,RApertures,EApertures,\
    irho,gap,A,B,L1,L2,K1,K2,max_order,num_int_steps,E0,scaling,\
    entrance_angle,useFringe1,fint1,h1,exit_angle,useFringe2,fint2,h2)
      for (int c = 0; c<num_particles; c++) { /* Loop over particles */
        double *r6 = r + 6*c;
        if (!atIsNaN(r6[0])) {
            int m;
            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, scaling);
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* edge focus */
            if (useFringe1)
                edge_fringe2A(r6, irho, entrance_angle,fint1,gap,h1,B[1]);
                /* edge_fringe(r6, irho, entrance_angle,fint1,gap); */
            else
                edge_fringe2A(r6, irho, entrance_angle,0,0,h1,B[1]);
            /* integrator */
            for (m=0; m < num_int_steps; m++) { /* Loop over slices */
                ATbendhxdrift6(r6,L1,irho);
                bndthinkickrad(r6, A, B, K1, irho, E0, max_order);
                ATbendhxdrift6(r6,L2,irho);
                bndthinkickrad(r6, A, B, K2, irho, E0, max_order);
                ATbendhxdrift6(r6,L2,irho);
                bndthinkickrad(r6, A, B, K1, irho, E0, max_order);
                ATbendhxdrift6(r6,L1,irho);
            }
            /* edge focus */
            if (useFringe2)
                edge_fringe2B(r6, irho, exit_angle, fint2, gap, h2, B[1]);
                /*	edge_fringe(r6, irho, exit_angle,fint2,gap);  */
            else
                edge_fringe2B(r6, irho, exit_angle, 0, 0, h2, B[1]);
            /* Check physical apertures at the exit of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);
            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, 1.0/scaling);
        }
    }
    if (KickAngle) {  /* Remove corrector component in polynomial coefficients */
        B[0] = B0;
        A[0] = A0;
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double irho, energy;
    if (!Elem) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, Scaling, Energy,
                FringeInt1, FringeInt2;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, h1, h2, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        Energy=atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        h1=atGetOptionalDouble(ElemData,"H1",0); check_error();
        h2=atGetOptionalDouble(ElemData,"H2",0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        Elem->NumIntSteps=NumIntSteps;
        Elem->BendingAngle=BendingAngle;
        Elem->EntranceAngle=EntranceAngle;
        Elem->ExitAngle=ExitAngle;
        Elem->Energy=Energy;
        /*optional fields*/
        Elem->FullGap=FullGap;
        Elem->Scaling=Scaling;
        Elem->FringeInt1=FringeInt1;
        Elem->FringeInt2=FringeInt2;
        Elem->h1=h1;
        Elem->h2=h2;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->KickAngle=KickAngle;
    }
    irho = Elem->BendingAngle/Elem->Length;
    energy = atEnergy(Param->energy, Elem->Energy);

    BndMPoleSymplectic4E2RadPass(r_in, Elem->Length, irho, Elem->PolynomA, Elem->PolynomB,
            Elem->MaxOrder, Elem->NumIntSteps, Elem->EntranceAngle, Elem->ExitAngle,
            Elem->FringeInt1, Elem->FringeInt2, Elem->FullGap,
            Elem->h1, Elem->h2,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            Elem->RApertures, Elem->EApertures,
            Elem->KickAngle, Elem->Scaling, energy, num_particles);
    return Elem;
}

MODULE_DEF(BndMPoleSymplectic4E2RadPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double rest_energy = 0.0;
        double charge = -1.0;
        double irho;
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, Scaling, FringeInt1, FringeInt2, Energy;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, h1, h2, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *KickAngle;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgTxt("Second argument must be a 6 x N matrix");
        
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        Energy=atGetOptionalDouble(ElemData,"Energy",0.0); check_error();
        FullGap=atGetOptionalDouble(ElemData,"FullGap", 0); check_error();
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1", 0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2", 0); check_error();
        h1=atGetOptionalDouble(ElemData,"H1", 0); check_error();
        h2=atGetOptionalDouble(ElemData,"H2", 0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();
        irho = BendingAngle/Length;
        if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        BndMPoleSymplectic4E2RadPass(r_in, Length, irho, PolynomA, PolynomB,
            MaxOrder, NumIntSteps, EntranceAngle, ExitAngle,
            FringeInt1, FringeInt2, FullGap,
            h1, h2,
            T1, T2, R1, R2, RApertures, EApertures,
            KickAngle, Scaling, Energy, num_particles);
    } else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(9,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("BendingAngle"));
        mxSetCell(plhs[0],2,mxCreateString("EntranceAngle"));
        mxSetCell(plhs[0],3,mxCreateString("ExitAngle"));
        mxSetCell(plhs[0],4,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],5,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],6,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],7,mxCreateString("NumIntSteps"));

        if (nlhs>1) {    /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(12,1);
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
            mxSetCell(plhs[1],11,mxCreateString("FieldScaling"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
