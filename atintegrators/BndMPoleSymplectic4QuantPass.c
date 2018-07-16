#include <math.h>
#include "atelem.c"
#include "atlalib.c"
#include "atphyslib.c"		/* edge, edge_fringe */
#include "driftkick.c"		/* fastdrift.c, bndthinkick.c */
#include "quadfringe.c"		/* QuadFringePassP, QuadFringePassN */
#include "atquantlib.c"

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

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
    int FringeBendEntrance;
    int FringeBendExit;
    double FringeInt1;
    double FringeInt2;
    double FullGap;
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
};

void BndMPoleSymplectic4QuantPass(double *r, double le, double irho, double *A, double *B,
        int max_order, int num_int_steps,
        double entrance_angle, double exit_angle,
        int FringeBendEntrance, int FringeBendExit,
        double fint1, double fint2, double gap,
        int FringeQuadEntrance, int FringeQuadExit,
        double *fringeIntM0,  /* I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 */
        double *fringeIntP0,  /* I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1 */
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double E0, int num_particles)
{
    double *r6;
    double dpp0, ng, ec, de, energy, gamma, cstec, cstng;
    double s0, ds, rho, dxp, dyp, xp0, yp0;
    int nph;
    int i, c, m;
    double p_norm, NormL1, NormL2;
    bool useT1 = (T1 != NULL);
    bool useT2 = (T2 != NULL);
    bool useR1 = (R1 != NULL);
    bool useR2 = (R2 != NULL);
    double SL = le/num_int_steps;
    double L1 = SL*DRIFT1;
    double L2 = SL*DRIFT2;
    double K1 = SL*KICK1;
    double K2 = SL*KICK2;
    bool useLinFrEleEntrance = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadEntrance==2);
    bool useLinFrEleExit = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadExit==2);
    double  qe = 1.60217733e-19;
    double  epsilon0 = 8.854187817e-12;
    double  clight = 2.99792458e8;
    double  emass = 510998.9461; /* electron mass in eV */ /* 9.10938188e-31; in kg*/
    double  hbar = 1.054571726e-34;
    double  pi = 3.14159265358979;
    double  alpha0 = qe*qe/(4*pi*epsilon0*hbar*clight);
    
    for (c = 0; c<num_particles; c++) {	/* Loop over particles  */
        r6 = r+c*6;
        if(!atIsNaN(r6[0])) {
            /*  misalignment at entrance  */
            if(useT1) ATaddvv(r6,T1);
            if(useR1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* edge focus */
            edge_fringe_entrance(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance);
            /* quadrupole gradient fringe */
            if (FringeQuadEntrance && B[1]!=0) {
                if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
                    linearQuadFringeElegantEntrance(r6, B[1], fringeIntM0, fringeIntP0);
                else
                    QuadFringePassP(r6, B[1]);
            }
            /* integrator */
            p_norm = 1/(1+r6[4]);
            NormL1 = L1*p_norm;
            NormL2 = L2*p_norm;
            for (m=0; m < num_int_steps; m++) {/* Loop over slices*/
                r6 = r+c*6;
                dpp0 = r6[4];
                xp0 = r6[1]/(1+r6[4]);
                yp0 = r6[3]/(1+r6[4]);
                s0 = r6[5];
                
                fastdrift(r6, NormL1);
                bndthinkick(r6, A, B, K1, irho, max_order);
                fastdrift(r6, NormL2);
                bndthinkick(r6, A, B, K2, irho, max_order);
                fastdrift(r6, NormL2);
                bndthinkick(r6, A, B,  K1, irho, max_order);
                fastdrift(r6, NormL1);
                
                energy = dpp0*E0+E0;
                
                gamma = energy/emass;/* emass in eV */
                cstec = 3.0*gamma*gamma*gamma*clight/(2.0)*hbar/qe;
                cstng = 5.0*sqrt(3.0)*alpha0*gamma/(6.0);
                
                dxp = r6[1]/(1+r6[4])-xp0 - irho*SL;
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
            /* quadrupole gradient fringe */
            if (FringeQuadExit && B[1]!=0) {
                if (useLinFrEleExit) /*Linear fringe fields from elegant*/
                    linearQuadFringeElegantExit(r6, B[1], fringeIntM0, fringeIntP0);
                else
                    QuadFringePassN(r6, B[1]);
            }
            /* edge focus */
            edge_fringe_exit(r6, irho, exit_angle, fint2, gap, FringeBendExit);
            
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
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap,
                FringeInt1, FringeInt2, Energy;
        int MaxOrder, NumIntSteps,  FringeBendEntrance, FringeBendExit,
                FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2,
                *EApertures, *RApertures, *fringeIntM0, *fringeIntP0;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
        FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        FringeQuadEntrance=atGetOptionalLong(ElemData,"FringeQuadEntrance",1); check_error();
        FringeQuadExit=atGetOptionalLong(ElemData,"FringeQuadExit",1); check_error();
        fringeIntM0=atGetOptionalDoubleArray(ElemData,"fringeIntM0"); check_error();
        fringeIntP0=atGetOptionalDoubleArray(ElemData,"fringeIntP0"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
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
        Elem->FringeInt1=FringeInt1;
        Elem->FringeInt2=FringeInt2;
        Elem->FringeBendEntrance=FringeBendEntrance;
        Elem->FringeBendExit=FringeBendExit;
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
    }
    irho = Elem->BendingAngle/Elem->Length;
    BndMPoleSymplectic4QuantPass(r_in,Elem->Length,irho,Elem->PolynomA,Elem->PolynomB,
            Elem->MaxOrder,Elem->NumIntSteps,Elem->EntranceAngle,Elem->ExitAngle,
            Elem->FringeBendEntrance,Elem->FringeBendExit,
            Elem->FringeInt1,Elem->FringeInt2,Elem->FullGap,
            Elem->FringeQuadEntrance,Elem->FringeQuadExit,
            Elem->fringeIntM0,Elem->fringeIntP0,Elem->T1,Elem->T2,
            Elem->R1,Elem->R2,Elem->RApertures,Elem->EApertures,Elem->Energy,num_particles);
    return Elem;
}

MODULE_DEF(BndMPoleSymplectic4QuantPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2 ) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap,
                FringeInt1, FringeInt2, Energy;
        int MaxOrder, NumIntSteps, FringeBendEntrance, FringeBendExit,
                FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0;
        double irho;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
        FringeBendEntrance=atGetOptionalLong(ElemData,"FringeBendEntrance",1); check_error();
        FringeBendExit=atGetOptionalLong(ElemData,"FringeBendExit",1); check_error();
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
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
        irho = BendingAngle/Length;
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        BndMPoleSymplectic4QuantPass(r_in, Length, irho, PolynomA, PolynomB,
                MaxOrder,NumIntSteps,EntranceAngle,ExitAngle,
                FringeBendEntrance,FringeBendExit,FringeInt1,FringeInt2,
                FullGap,FringeQuadEntrance,FringeQuadExit,fringeIntM0,fringeIntP0,
                T1,T2,R1,R2,RApertures,EApertures,Energy,num_particles);
    }
    else if (nrhs == 0) {
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
        mxSetCell(plhs[0],8,mxCreateString("Energy"));
        
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(15,1);
            mxSetCell(plhs[1],0,mxCreateString("FullGap"));
            mxSetCell(plhs[1],1,mxCreateString("FringeInt1"));
            mxSetCell(plhs[1],2,mxCreateString("FringeInt2"));
            mxSetCell(plhs[1],3,mxCreateString("FringeBendEntrance"));
            mxSetCell(plhs[1],4,mxCreateString("FringeBendExit"));
            mxSetCell(plhs[1],5,mxCreateString("FringeQuadEntrance"));
            mxSetCell(plhs[1],6,mxCreateString("FringeQuadExit"));
            mxSetCell(plhs[1],7,mxCreateString("fringeIntM0"));
            mxSetCell(plhs[1],8,mxCreateString("fringeIntP0"));
            mxSetCell(plhs[1],9,mxCreateString("T1"));
            mxSetCell(plhs[1],10,mxCreateString("T2"));
            mxSetCell(plhs[1],11,mxCreateString("R1"));
            mxSetCell(plhs[1],12,mxCreateString("R2"));
            mxSetCell(plhs[1],13,mxCreateString("RApertures"));
            mxSetCell(plhs[1],14,mxCreateString("EApertures"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
