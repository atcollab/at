#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "atphyslib.c"
#include "driftkick.c"  /* strthinkick.c */

/* Straight dipole w/ multipole using Symplectic Integration and rotation at
 * dipole faces.
 * Created by Xiaobiao Huang, 7/31/2018 */
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
    /* Optional fields */
    double FringeInt1;
    double FringeInt2;
    double FullGap;
    double Scaling;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *KickAngle;
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
    /*double fx = inv_rho*tan(edge_angle);*/
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

void BndStrMPoleSymplectic4Pass(double *r, double le, double irho, double *A, double *B,
        int max_order, int num_int_steps,
        double entrance_angle, double exit_angle,
        double X0ref, double ByError, double RefDZ,
        double fint1, double fint2, double gap,
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double *KickAngle, double scaling, int num_particles)
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
    B[0] += irho;

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
	            edgey_fringe(r6, irho+B[1]*X0ref, entrance_angle,fint1,gap);
	        else
	            edgey(r6, irho+B[1]*X0ref, entrance_angle);
            /* Rotate and translate to straight Cartesian coordinate */
            E1rotation(r6, X0ref, entrance_angle);
            /* integrator */
            for (m=0; m < num_int_steps; m++) { /* Loop over slices */
				ladrift6(r6,L1);
			    strthinkick(r6, A, B, K1, max_order);
				ladrift6(r6,L2);
			    strthinkick(r6, A, B, K2, max_order);
				ladrift6(r6,L2);
				strthinkick(r6, A, B, K1, max_order);
				ladrift6(r6,L1);
			}
            /* Rotate and translate back to curvilinear coordinate */
            E2rotation(r6, X0ref, exit_angle);
            r6[5] -= RefDZ;
            /* edge focus */
			if (useFringe2)
	            edgey_fringe(r6, irho+B[1]*X0ref, exit_angle,fint2,gap);
	        else    /* edge focus */
	            edgey(r6, irho+B[1]*X0ref, exit_angle);
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
    B[0] = B0;
    A[0] = A0;
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double irho, flen;
    if (!Elem) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, Scaling,
                FringeInt1, FringeInt2, X0ref, ByError, RefDZ;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();
        X0ref=atGetOptionalDouble(ElemData,"X0ref",0); check_error();
        ByError=atGetOptionalDouble(ElemData,"ByError",0); check_error();
        RefDZ=atGetOptionalDouble(ElemData,"RefDZ",0); check_error();

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        Elem->NumIntSteps=NumIntSteps;
        Elem->BendingAngle=BendingAngle;
        Elem->EntranceAngle=EntranceAngle;
        Elem->ExitAngle=ExitAngle;
        /*optional fields*/
        Elem->FullGap=FullGap;
        Elem->Scaling=Scaling;
        Elem->FringeInt1=FringeInt1;
        Elem->FringeInt2=FringeInt2;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->KickAngle=KickAngle;
        Elem->X0ref=X0ref;
        Elem->ByError=ByError;
        Elem->RefDZ=RefDZ;
    }
    irho = Elem->BendingAngle / Elem->Length;
    flen = 2.0 / irho * sin(Elem->BendingAngle/2.0);
    BndStrMPoleSymplectic4Pass(r_in, flen, irho, Elem->PolynomA, Elem->PolynomB,
            Elem->MaxOrder, Elem->NumIntSteps, Elem->EntranceAngle, Elem->ExitAngle,
            Elem->X0ref, Elem->ByError, Elem->RefDZ,
            Elem->FringeInt1, Elem->FringeInt2, Elem->FullGap,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            Elem->RApertures, Elem->EApertures,
            Elem->KickAngle, Elem->Scaling, num_particles);
    return Elem;
}

MODULE_DEF(BndStrMPoleSymplectic4Pass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, FullGap, Scaling,
                FringeInt1, FringeInt2, X0ref, ByError, RefDZ;
        int MaxOrder, NumIntSteps;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures, *KickAngle;
        double irho, flen;
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
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        X0ref=atGetOptionalDouble(ElemData,"X0ref", 0); check_error();
        ByError=atGetOptionalDouble(ElemData,"ByError", 0); check_error();
        RefDZ=atGetOptionalDouble(ElemData,"RefDZ", 0); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();
        irho = BendingAngle/Length;
        flen = 2.0/irho*sin(BendingAngle/2.0); /* field length */

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        BndStrMPoleSymplectic4Pass(r_in, flen, irho, PolynomA, PolynomB,
            MaxOrder, NumIntSteps, EntranceAngle, ExitAngle,
            X0ref, ByError, RefDZ,
            FringeInt1, FringeInt2, FullGap,
            T1, T2, R1, R2, RApertures, EApertures,
            KickAngle, Scaling, num_particles);
    } else if (nrhs == 0) {
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

        if (nlhs>1) {    /* list of optional fields */
	        plhs[1] = mxCreateCellMatrix(14,1);
            mxSetCell(plhs[1],0,mxCreateString("FullGap"));
            mxSetCell(plhs[1],1,mxCreateString("FringeInt1"));
            mxSetCell(plhs[1],2,mxCreateString("FringeInt2"));
            mxSetCell(plhs[1],3,mxCreateString("X0ref"));
            mxSetCell(plhs[1],4,mxCreateString("ByError"));
            mxSetCell(plhs[1],5,mxCreateString("RefDZ"));
            mxSetCell(plhs[1],6,mxCreateString("T1"));
            mxSetCell(plhs[1],7,mxCreateString("T2"));
            mxSetCell(plhs[1],8,mxCreateString("R1"));
            mxSetCell(plhs[1],9,mxCreateString("R2"));
            mxSetCell(plhs[1],10,mxCreateString("RApertures"));
            mxSetCell(plhs[1],11,mxCreateString("EApertures"));
            mxSetCell(plhs[1],12,mxCreateString("KickAngle"));
            mxSetCell(plhs[1],13,mxCreateString("FieldScaling"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
