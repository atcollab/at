#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "atquantlib.c"
#include "driftkick.c"  	/* fastdrift.c, strthinkick.c */
#include "quadfringe.c"		/* QuadFringePassP, QuadFringePassN */

struct elem
{
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    double Energy;
    /* Optional fields */
    double Scaling;
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
        double *KickAngle, double scaling, double E0,
        pcg32_random_t* rng, int num_particles)
{
    double SL = le/num_int_steps;
    double L1 = SL*DRIFT1;
    double L2 = SL*DRIFT2;
    double K1 = SL*KICK1;
    double K2 = SL*KICK2;
    bool useLinFrEleEntrance = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadEntrance==2);
    bool useLinFrEleExit = (fringeIntM0 != NULL && fringeIntP0 != NULL  && FringeQuadExit==2);
    double qe = 1.60217733e-19;
    double epsilon0 = 8.854187817e-12;
    double clight = 2.99792458e8;
    double emass = 510998.9461; /* electron mass in eV */ /* 9.10938188e-31; in kg*/
    double hbar = 1.054571726e-34;
    double pi = 3.14159265358979;
    double alpha0 = qe * qe / (4 * pi * epsilon0 * hbar * clight);
    double B0 = B[0];
    double A0 = A[0];

    if (KickAngle) {   /* Convert corrector component to polynomial coefficients */
        B[0] -= sin(KickAngle[0])/le;
        A[0] += sin(KickAngle[1])/le;
    }
/*  The behaviour of random generators with OpenMP is doubtful. OpenMP disabled until
    it's understood
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none)              \
    shared(r,num_particles,R1,T1,R2,T2,RApertures,EApertures,                                       \
    A,B,L1,L2,K1,K2,max_order,num_int_steps,rng,scaling,                                            \
    FringeQuadEntrance, useLinFrEleEntrance,FringeQuadExit,useLinFrEleExit,fringeIntM0,fringeIntP0, \
    emass,E0,hbar,clight,alpha0,qe,SL)
*/
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
            if (FringeQuadEntrance && B[1]!=0) {
                if (useLinFrEleEntrance) /*Linear fringe fields from elegant*/
                    linearQuadFringeElegantEntrance(r6, B[1], fringeIntM0, fringeIntP0);
                else
                    QuadFringePassP(r6, B[1]);
            }
            /* integrator */
            for (m=0; m < num_int_steps; m++) { /* Loop over slices */
                int i;
                double ng, ec, de, energy, gamma, cstec, cstng;
                double ds, rho, dxp, dyp;
                int nph;
                double p_norm = 1.0 / (1.0 + r6[4]);
                double NormL1 = L1 * p_norm;
                double NormL2 = L2 * p_norm;
                double dpp0 = r6[4];
                double xp0 = r6[1] * p_norm;
                double yp0 = r6[3] * p_norm;
                double s0 = r6[5];

                fastdrift(r6, NormL1);
                strthinkick(r6, A, B, K1, max_order);
                fastdrift(r6, NormL2);
                strthinkick(r6, A, B, K2, max_order);
                fastdrift(r6, NormL2);
                strthinkick(r6, A, B, K1, max_order);
                fastdrift(r6, NormL1);

                energy = dpp0 * E0 + E0;

                gamma = energy / emass; /* emass in eV */
                cstec = 3.0 * gamma * gamma * gamma * clight / (2.0) * hbar / qe;
                cstng = 5.0 * sqrt(3.0) * alpha0 * gamma / (6.0);

                dxp = r6[1] * p_norm - xp0;
                dyp = r6[3] * p_norm - yp0;
                ds = r6[5] - s0;

                rho = (SL + ds) / sqrt(dxp * dxp + dyp * dyp);

                ng = cstng / rho * (SL + ds);
                ec = cstec / rho;

                nph = atrandp_r(rng, ng);

                de = 0.0;
                for (i = 0; i < nph; i++) {
                    de = de + getEnergy(rng, ec);
                };
                r6[4] = r6[4] - de / E0;
                r6[1] = r6[1] * p_norm * (1 + r6[4]);
                r6[3] = r6[3] * p_norm * (1 + r6[4]);
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
    double energy;
    if (!Elem) {
        double Length, Energy, Scaling;
        int MaxOrder, NumIntSteps, FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0, *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        /*optional fields*/
        Energy=atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
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

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        Elem->NumIntSteps=NumIntSteps;
        Elem->Energy=Energy;
        /*optional fields*/
        Elem->Scaling=Scaling;
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
    energy = atEnergy(Param->energy, Elem->Energy);

    StrMPoleSymplectic4QuantPass(r_in, Elem->Length, Elem->PolynomA, Elem->PolynomB,
            Elem->MaxOrder, Elem->NumIntSteps,
            Elem->FringeQuadEntrance, Elem->FringeQuadExit,
            Elem->fringeIntM0, Elem->fringeIntP0,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            Elem->RApertures, Elem->EApertures,
            Elem->KickAngle, Elem->Scaling,
            energy, Param->thread_rng, num_particles);
    return Elem;
}

MODULE_DEF(StrMPoleSymplectic4QuantPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double rest_energy = 0.0;
        double charge = -1.0;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length, Energy, Scaling;
        int MaxOrder, NumIntSteps, FringeQuadEntrance, FringeQuadExit;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *fringeIntM0, *fringeIntP0, *KickAngle;
        if (mxGetM(prhs[1]) != 6) mexErrMsgTxt("Second argument must be a 6 x N matrix");

        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        /*optional fields*/
        Energy=atGetOptionalDouble(ElemData,"Energy",0.0); check_error();
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
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
        if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);

        StrMPoleSymplectic4QuantPass(r_in, Length, PolynomA, PolynomB,
            MaxOrder, NumIntSteps,
            FringeQuadEntrance, FringeQuadExit,
            fringeIntM0, fringeIntP0,
            T1, T2, R1, R2, RApertures, EApertures,
            KickAngle, Scaling, Energy, &pcg32_global, num_particles);
    } else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(6, 1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],2,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],3,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],4,mxCreateString("NumIntSteps"));
        mxSetCell(plhs[0],5,mxCreateString("Energy"));
        if (nlhs>1) {    /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(12,1);
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
            mxSetCell(plhs[1],11,mxCreateString("FieldScaling"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
