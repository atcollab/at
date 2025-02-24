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
    double HarmNumber;
    double TimeLag;
    double PhaseLag;
    double SigPhi;
    double SigVV;
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

void CrabCavityPass(
    double *r_in,
    double le,
    double nvx,
    double nvy,
    double freq,
    double h,
    double lag,
    double philag,
    int nturn,
    double T0,
    double gamma0,
    int num_particles)
{
    double beta0, betgam0, phiref;
    double k = TWOPI*freq/C0;
    kickfun ff;

    if (isfinite(gamma0)) {
        betgam0 = sqrt(gamma0*gamma0 - 1.0);
        beta0 = betgam0 / gamma0;
        ff = &non_relativistic_kick;
        nvx = nvx/beta0/beta0;
        nvy = nvy/beta0/beta0;
    }
    else {
        betgam0 = gamma0;
        beta0 = 1.0;
        ff = &relativistic_kick;
    }
    phiref = TWOPI*(freq*T0 - h)*nturn - k*lag/beta0 - philag;

    if (le == 0.0) {
        #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) \
            default(none) shared(r_in, num_particles, normv, k, betgam0, phiref)
        for (int c = 0; c < num_particles; c++) {
            double *r6 = r_in+c*6;
            if (!atIsNaN(r6[0]))
                /* Longitudinal momentum kick */
                ff(r6, nvx, nvy, k, betgam0, phiref);
        }
    }
    else {
        double halflength = le/2;
        #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) \
            default(none) shared(r_in, num_particles, normv, k, betgam0, phiref, halflength)
        for (int c = 0; c < num_particles; c++) {
            double *r6 = r_in+c*6;
            if (!atIsNaN(r6[0]))  {
                /* Propagate through a drift equal to half cavity length */
                drift6(r6, halflength);
                /* Longitudinal momentum kick */
                ff(r6, nvx, nvy, k, betgam0, phiref);
                /* Propagate through a drift equal to half cavity length */
                drift6(r6, halflength);
            }
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{
    double energy, nvx, nvy;
    double PhaseLag, SigPhi, SigVV;

    if (!Elem) {
        double Length, Frequency, *Voltages, Energy, HarmNumber, TimeLag;
		Length = atGetDouble(ElemData,"Length"); check_error();
		Voltages = atGetDoubleArray(ElemData,"Voltages"); check_error();
		Frequency = atGetDouble(ElemData,"Frequency"); check_error();
		/* Optional fields */
		Energy = atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
        HarmNumber=atGetOptionalLong(ElemData,"HarmNumber", round(Frequency*Param->T0)); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0.0); check_error();
        PhaseLag=atGetOptionalDouble(ElemData,"PhaseLag",0.0); check_error();
		SigPhi = atGetOptionalDouble(ElemData,"SigPhi",0.0); check_error();
		SigVV = atGetOptionalDouble(ElemData,"SigVV",0.0); check_error();

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->Vx=Voltages[0];
        Elem->Vy=Voltages[1];
        Elem->Energy=Energy;
        Elem->Frequency=Frequency;
        Elem->HarmNumber=HarmNumber;
        Elem->TimeLag=TimeLag;
        Elem->PhaseLag=PhaseLag;
        Elem->SigPhi=SigPhi;
        Elem->SigVV=SigVV;
    }
    else {
        PhaseLag=Elem->PhaseLag;
        SigPhi=Elem->SigPhi;
        SigVV=Elem->SigVV;
    }
    energy = atEnergy(Param->energy, Elem->Energy);
    nvx = Elem->Vx;
    nvy = Elem->Vy;

    if (SigPhi > 0.0)
        PhaseLag += atrandn_r(Param->common_rng, 0.0, SigPhi);
    if (SigVV > 0.0) {
        nvx *= atrandn_r(Param->common_rng, 1.0, SigVV);
        nvy *= atrandn_r(Param->common_rng, 1.0, SigVV);
    }

    CrabCavityPass(r_in, Elem->Length, nvx/energy, nvy/energy, Elem->Frequency, Elem->HarmNumber,
        Elem->TimeLag, PhaseLag, Param->nturn, Param->T0, energy/Param->rest_energy, num_particles);
    return Elem;
}

MODULE_DEF(CrabCavityPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double nvx, nvy;
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
        double HarmNumber=atGetOptionalLong(ElemData,"HarmNumber", 0); check_error();
		double TimeLag = atGetOptionalDouble(ElemData,"TimeLag",0.0); check_error();
		double PhaseLag = atGetOptionalDouble(ElemData,"PhaseLag",0.0); check_error();
        if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);

        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        nvx = Voltages[0] / Energy;
        nvy = Voltages[1] / Energy;
	    CrabCavityPass(r_in, Length, nvx, nvy, Frequency, HarmNumber, TimeLag, PhaseLag,
	        0, 0.0, Energy/rest_energy, num_particles);
    }
    else if (nrhs == 0) {
	    plhs[0] = mxCreateCellMatrix(3,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("Voltages"));
	    mxSetCell(plhs[0],2,mxCreateString("Frequency"));
	    if (nlhs > 1) { /* optional fields */
	        plhs[1] = mxCreateCellMatrix(6,1);
    	    mxSetCell(plhs[0],0,mxCreateString("Energy"));
    	    mxSetCell(plhs[0],1,mxCreateString("HarmNumber"));
            mxSetCell(plhs[1],2,mxCreateString("TimeLag"));
            mxSetCell(plhs[1],3,mxCreateString("PhaseLag"));
			mxSetCell(plhs[1],4,mxCreateString("SigPhi"));
			mxSetCell(plhs[1],5,mxCreateString("SigVV"));
	    }
    }
    else {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
