/* CrabCavityPass.c
   for Accelerator Toolbox 
   Xiaobiao Huang, created 11/1/2016
*/

#include "atelem.c"
#include "atrandom.c"
#include "driftkickrad.c"

struct elem
{
    double Length;
    double Vx;
    double Vy;
    double Frequency;
    double Energy;
    double TimeLag;
    double PhaseLag;
    double SigPhi;
    double SigVV;
    long HarmNumber;
};

typedef void (*kickfun)(double *r_in, double nvx, double nvy, double k, double betgam0, double phi);

static void relativistic_kick(double *r_in, double nvx, double nvy, double k, double betgam0, double phi)
{
    double omega_t = k*r_in[5] + phi;
    double cos_ot = cos(omega_t);
    double sin_ot = sin(omega_t);
    double ddpp = -(nvy*r_in[2]+nvx*r_in[0])*k*cos_ot/(1.0+r_in[4]);

    r_in[1] += nvx*sin_ot;
    r_in[3] += nvy*sin_ot;
    r_in[4] += ddpp;
}

static void non_relativistic_kick(double *r_in, double nvx, double nvy, double k, double betgam0, double phi)
{
    double betgami = betgam0*(1.0 + r_in[4]);
    double betai = betgami/sqrt(1.0 + betgami*betgami);
    double omega_t = k*r_in[5]/betai + phi;
    double cos_ot = cos(omega_t);
    double sin_ot = sin(omega_t);
    double ddpp = -(nvy*r_in[2]+nvx*r_in[0])*k*cos_ot/(1.0+r_in[4]);

    r_in[1] += nvx*sin_ot;
    r_in[3] += nvy*sin_ot;
    r_in[4] += ddpp;
}

static void CrabCavityPass(
    double *r_in,
    double le,
    double nvx,
    double nvy,
    double freq,
    double betgam0,
    double phiref,
    int num_particles)
{
	double k = TWOPI*freq/C0;
    kickfun ff = isfinite(betgam0) ? &non_relativistic_kick : &relativistic_kick;

	if (le == 0) {
        #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) \
            default(none) shared(r_in, num_particles, nvx, nvy, k, phiref)
        for (int c = 0; c<num_particles; c++) {
            double *r6 = r_in + c*6;
            if (!atIsNaN(r6[0])) {
                ff(r6, nvx, nvy, k, betgam0, phiref);
            }
        }
	}
	else {
	    double halflength = le/2.0;

        #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) \
            default(none) shared(r_in, num_particles, nvx, nvy, k, phiref, halflength)
        for (int c = 0; c<num_particles; c++) {
            double *r6 = r_in + c*6;
            if (!atIsNaN(r6[0])) {
                drift6(r6, halflength);
                ff(r6, nvx, nvy, k, betgam0, phiref);
                drift6(r6, halflength);
            }
        }
    }

}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{
    double energy, beta0, gamma0, betgam0;
    double t0f, lag, phiref, nvx, nvy;
    long HarmNumber;
    double SigPhi, SigVV;

    if (!Elem) {
        double Length, *Voltages, Energy, Frequency, TimeLag, PhaseLag;
		Length = atGetDouble(ElemData,"Length"); check_error();
		Voltages = atGetDoubleArray(ElemData,"Voltages"); check_error();
		Frequency = atGetDouble(ElemData,"Frequency"); check_error();
		/* Optional fields */
		Energy = atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0.0); check_error();
        PhaseLag=atGetOptionalDouble(ElemData,"PhaseLag",0.0); check_error();
        HarmNumber=atGetOptionalLong(ElemData,"HarmNumber", lround(Frequency*Param->T0)); check_error();
		SigPhi = atGetOptionalDouble(ElemData,"SigPhi",0.0); check_error();
		SigVV = atGetOptionalDouble(ElemData,"SigVV",0.0); check_error();

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->Vx=Voltages[0];
        Elem->Vy=Voltages[1];
        Elem->Energy=Energy;
        Elem->Frequency=Frequency;
        Elem->TimeLag=TimeLag;
        Elem->PhaseLag=PhaseLag;
        Elem->SigPhi=SigPhi;
        Elem->SigVV=SigVV;
    }
    else {
        SigPhi=Elem->SigPhi;
        SigVV=Elem->SigVV;
    }
    energy = atEnergy(Param->energy, Elem->Energy);
    gamma0 = energy/Param->rest_energy;

    if (isfinite(gamma0)) {
        betgam0 = sqrt(gamma0*gamma0 - 1.0);
        beta0 = betgam0 / gamma0;
    }
    else {
        betgam0 = gamma0;
        beta0 = 1.0;
    }
    t0f = Elem->Frequency * Param->T0;
    lag = TWOPI*Elem->Frequency*Elem->TimeLag/beta0/C0 + Elem->PhaseLag;
    phiref = TWOPI * (t0f - Elem->HarmNumber) * Param->nturn - lag;
    nvx = Elem->Vx/energy/beta0/beta0;
    nvy = Elem->Vy/energy/beta0/beta0;

    if (SigPhi > 0.0)
        phiref += atrandn_r(Param->common_rng, 0.0, SigPhi);
    if (SigVV > 0.0) {
        nvx *= atrandn_r(Param->common_rng, 1.0, SigVV);
        nvy *= atrandn_r(Param->common_rng, 1.0, SigVV);
    }

    CrabCavityPass(r_in, Elem->Length, nvx, nvy, Elem->Frequency, betgam0, phiref, num_particles);
    return Elem;
}

MODULE_DEF(CrabCavityPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double betgam0, beta0, gamma0;
        double lag, phiref;
        double *r_in;
        double rest_energy = 0.0;
        double charge = -1.0;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length=atGetDouble(ElemData,"Length");
		double *Voltages = atGetDoubleArray(ElemData,"Voltages"); check_error();
		double Frequency = atGetDouble(ElemData,"Frequency"); check_error();
		/* Optional fields */
		double Energy = atGetOptionalDouble(ElemData,"Energy",0.0); check_error();
		double TimeLag = atGetOptionalDouble(ElemData,"TimeLag",0.0); check_error();
		double PhaseLag = atGetOptionalDouble(ElemData,"PhaseLag",0.0); check_error();
        if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);
        gamma0 = Energy/rest_energy;

        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        if (isfinite(gamma0)) {
            betgam0 = sqrt(gamma0*gamma0 - 1.0);
            beta0 = betgam0 / gamma0;
        }
        else {
            betgam0 = gamma0;
            beta0 = 1.0;
        }
        lag = TWOPI*Frequency*TimeLag/beta0/C0 + PhaseLag;
        phiref = -lag;
        normvx = Voltages[0] / Energy / beta0 / beta0;
        normvy = Voltages[1] / Energy / beta0 / beta0;
	    CrabCavityPass(r_in, Length, normvx, normvy, Frequency, betgam0, phiref, num_particles);
    }
    else if (nrhs == 0) {
	    plhs[0] = mxCreateCellMatrix(3,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("Voltages"));
	    mxSetCell(plhs[0],2,mxCreateString("Frequency"));
	    if (nlhs > 1) { /* optional fields */
	        plhs[1] = mxCreateCellMatrix(5,1);
    	    mxSetCell(plhs[0],0,mxCreateString("Energy"));
            mxSetCell(plhs[1],1,mxCreateString("TimeLag"));
            mxSetCell(plhs[1],2,mxCreateString("PhaseLag"));
			mxSetCell(plhs[1],3,mxCreateString("SigPhi"));
			mxSetCell(plhs[1],4,mxCreateString("SigVV"));
	    }
    }
    else {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
