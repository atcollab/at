#include "atelem.c"
#include "atlalib.c"
#include "driftkickrad.c"	/* drift6.c, strthinkickrad.c */
#include "quadfringe.c"		/* QuadFringePassP, QuadFringePassN */
#include "atquantlib.c"

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

#define SQR(X) ((X)*(X))

struct elem
{
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    double Energy;
    /* Optional fields */
    int FringeQuadEntrance;
    int FringeQuadExit;
    double *fringeIntM0;
    double *fringeIntP0;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *KickAngle;
};

void StrMPoleSymplectic4QuantPass(double *r, double le, double *A, double *B,
        int max_order, int num_int_steps,
        int FringeQuadEntrance, int FringeQuadExit, /* 0 (no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) */
        double *fringeIntM0,  /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
        double *fringeIntP0,  /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double *KickAngle, double E0,
        int num_particles)
{
    double *r6;
    double SL, L1, L2, K1, K2;
    double dpp0, ng, ec, de, energy, gamma, cstec, cstng;
    double s0, ds, rho, dxp, dyp, xp0, yp0;
    int c,m,i;
    int nph;
    double  qe = 1.60217733e-19;
    double  epsilon0 = 8.854187817e-12;
    double  clight = 2.99792458e8;
    double  emass = 510998.9461; /* electron mass in eV */ /* 9.10938188e-31; in kg*/
    double  hbar = 1.054571726e-34;
    double  pi = 3.14159265358979;
    double  alpha0 = qe*qe/(4*pi*epsilon0*hbar*clight);
    bool useLinFrEleEntrance = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadEntrance==2);
    bool useLinFrEleExit = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadExit==2);
    SL = le/num_int_steps;
    L1 = SL*DRIFT1;
    L2 = SL*DRIFT2;
    K1 = SL*KICK1;
    K2 = SL*KICK2;
    
    if (KickAngle) {  /* Convert corrector component to polynomial coefficients */
        B[0] -= sin(KickAngle[0])/le;
        A[0] += sin(KickAngle[1])/le;
    }
#pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(shared) shared(r,num_particles) private(c,r6,m)
    for (c = 0;c<num_particles;c++)	{   /* Loop over particles  */
        r6 = r+c*6;
        if(!atIsNaN(r6[0])) {
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            if (FringeQuadEntrance && B[1]!=0) {
                if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
                    linearQuadFringeElegantEntrance(r6, B[1], fringeIntM0, fringeIntP0);
                else
                    QuadFringePassP(r6, B[1]);
            }
            /* integrator */
            for (m=0; m < num_int_steps; m++) { /* Loop over slices */
                r6 = r+c*6;
                dpp0 = r6[4];
                xp0 = r6[1]/(1+r6[4]);
                yp0 = r6[3]/(1+r6[4]);
                s0 = r6[5];
                
                ATdrift6(r6,L1);
                strthinkickrad(r6, A, B, K1, E0, max_order);
                ATdrift6(r6,L2);
                strthinkickrad(r6, A, B, K2, E0, max_order);
                ATdrift6(r6,L2);
                strthinkickrad(r6, A, B, K1, E0, max_order);
                ATdrift6(r6,L1);
                
                energy = dpp0*E0+E0;
                
                gamma = energy/emass;/* emass in eV */
                cstec = 3.0*gamma*gamma*gamma*clight/(2.0)*hbar/qe;
                cstng = 5.0*sqrt(3.0)*alpha0*gamma/(6.0);
                
                dxp = r6[1]/(1+r6[4])-xp0;
                dyp = r6[3]/(1+r6[4])-yp0;
                ds = r6[5]-s0;
                
                rho = (SL+ds)/sqrt(dxp*dxp+dyp*dyp);
                
                ng =  cstng*1/rho*(SL+ds);
                ec =  cstec*1/rho;
                
                nph = poissonRandomNumber(ng);
                de = 0.0;
                if(nph>0){
                    for(i=0;i<nph;i++){
                        de = de + getEnergy(ec);
                    };
                };
                
                r6[1] = r6[1]/(1+r6[4])*(1+r6[4]-de/E0);
                r6[3] = r6[3]/(1+r6[4])*(1+r6[4]-de/E0);
                r6[4] = r6[4]-de/E0;
            }
            if (FringeQuadExit && B[1]!=0) {
                if (useLinFrEleExit) /*Linear fringe fields from elegant*/
                    linearQuadFringeElegantExit(r6, B[1], fringeIntM0, fringeIntP0);
                else
                    QuadFringePassN(r6, B[1]);
            }
            /* Check physical apertures at the exit of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);
        }
    }
    if (KickAngle) {  /* Remove corrector component in polynomial coefficients */
        B[0] += sin(KickAngle[0])/le;
        A[0] -= sin(KickAngle[1])/le;
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Length, Energy;
        int MaxOrder, NumIntSteps, FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0, *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
        FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0);
        FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0);
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
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
        Elem->Energy=Energy;
        /*optional fields*/
        Elem->FringeQuadEntrance=FringeQuadEntrance;
        Elem->FringeQuadExit=FringeQuadExit;
        Elem->fringeIntM0=fringeIntM0;
        Elem->fringeIntP0=fringeIntP0;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->KickAngle=KickAngle;
    }
    StrMPoleSymplectic4QuantPass(r_in,Elem->Length,Elem->PolynomA,Elem->PolynomB,
            Elem->MaxOrder,Elem->NumIntSteps,Elem->FringeQuadEntrance,
            Elem->FringeQuadExit,Elem->fringeIntM0,Elem->fringeIntP0,
            Elem->T1,Elem->T2,Elem->R1,Elem->R2,
            Elem->RApertures,Elem->EApertures,Elem->KickAngle,Elem->Energy,num_particles);
    return Elem;
}

MODULE_DEF(StrMPoleSymplectic4QuantPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length, Energy;
        int MaxOrder, NumIntSteps, FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0, *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
        FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",0); check_error();
        FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",0); check_error();
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        StrMPoleSymplectic4QuantPass(r_in,Length,PolynomA,PolynomB,MaxOrder,NumIntSteps,
                FringeQuadEntrance,FringeQuadExit,fringeIntM0,fringeIntP0,
                T1,T2,R1,R2,RApertures,EApertures,KickAngle,Energy,num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],2,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],3,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],4,mxCreateString("NumIntSteps"));
        mxSetCell(plhs[0],5,mxCreateString("Energy"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(11,1);
            mxSetCell(plhs[1], 0,mxCreateString("FringeQuadEntrance"));
            mxSetCell(plhs[1], 1,mxCreateString("FringeQuadExit"));
            mxSetCell(plhs[1], 2,mxCreateString("fringeIntM0"));
            mxSetCell(plhs[1], 3,mxCreateString("fringeIntP0"));
            mxSetCell(plhs[1], 4,mxCreateString("T1"));
            mxSetCell(plhs[1], 5,mxCreateString("T2"));
            mxSetCell(plhs[1], 6,mxCreateString("R1"));
            mxSetCell(plhs[1], 7,mxCreateString("R2"));
            mxSetCell(plhs[1], 8,mxCreateString("RApertures"));
            mxSetCell(plhs[1], 9,mxCreateString("EApertures"));
            mxSetCell(plhs[1],10,mxCreateString("KickAngle"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/

