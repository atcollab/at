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
    double Phase;
    double SigPhi;
    double SigVV;
};

static void track_cavity(double *r_in, double nvx, double nvy, double k, double phi)
{
    double omega_t = k*r_in[5] + phi;
    double cos_ot = cos(omega_t);
    double sin_ot = sin(omega_t);
    double ddpp = -(nvy*r_in[2]+nvx*r_in[0])*k*cos_ot/(1.0+r_in[4]);

    r_in[1] += nvx*sin_ot;
    r_in[3] += nvy*sin_ot;
    r_in[4] += ddpp;
}

static void CrabCavityPass(double *r_in, double le, double nvx, double nvy,
    double freq, double phi, double sigphi, double sigvv, pcg32_random_t* rng, int num_particles)
{
	double k = TWOPI*freq/C0;

    if (sigphi > 0.0)
        phi = phi + atrandn_r(rng, 0.0, sigphi);
    if (sigvv > 0.0) {
        nvx *= atrandn_r(rng, 1.0, sigvv);
        nvy *= atrandn_r(rng, 1.0, sigvv);
    }

	if (le == 0) {
        #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) \
            default(none) shared(r_in, num_particles, nvx, nvy, k, phi)
        for (int c = 0; c<num_particles; c++) {
            double *r6 = r_in + c*6;
            if (!atIsNaN(r6[0])) {
                track_cavity(r6, nvx, nvy, k, phi);
            }
        }
	}
	else {
	    double halflength = le/2.0;

        #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) \
            default(none) shared(r_in, num_particles, nvx, nvy, k, phi, halflength)
        for (int c = 0; c<num_particles; c++) {
            double *r6 = r_in + c*6;
            if (!atIsNaN(r6[0])) {
                drift6(r6, halflength);
                track_cavity(r6, nvx, nvy, k, phi);
                drift6(r6, halflength);
            }
        }
    }

}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{
    double energy;

    if (!Elem) {
        double Length, *Voltages, Energy, Frequency, Phase, SigPhi, SigVV;
		Length = atGetDouble(ElemData,"Length"); check_error();
		Voltages = atGetDoubleArray(ElemData,"Voltages"); check_error();
		Frequency = atGetDouble(ElemData,"Frequency"); check_error();
		/* Optional fields */
		Energy = atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
		Phase = atGetOptionalDouble(ElemData,"Phase",0.0); check_error();
		SigPhi = atGetOptionalDouble(ElemData,"SigPhi",0.0); check_error();
		SigVV = atGetOptionalDouble(ElemData,"SigVV",0.0); check_error();

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->Vx=Voltages[0];
        Elem->Vy=Voltages[1];
        Elem->Energy=Energy;
        Elem->Frequency=Frequency;
        Elem->Phase=Phase;
        Elem->SigPhi=SigPhi;
        Elem->SigVV=SigVV;
    }
    energy = atEnergy(Param->energy, Elem->Energy);
    CrabCavityPass(r_in, Elem->Length, Elem->Vx/energy, Elem->Vy/energy,
        Elem->Frequency, Elem->Phase, Elem->SigPhi, Elem->SigVV, Param->common_rng, num_particles);
    return Elem;
}

MODULE_DEF(CrabCavityPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
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
		double Phase = atGetOptionalDouble(ElemData,"Phase",0.0); check_error();
		double SigPhi = atGetOptionalDouble(ElemData,"SigPhi",0.0); check_error();
		double SigVV = atGetOptionalDouble(ElemData,"SigVV",0.0); check_error();
        if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);

        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
	    CrabCavityPass(r_in, Length, Voltages[0]/Energy, Voltages[1]/Energy,
	        Frequency, Phase, SigPhi, SigVV, &pcg32_global, num_particles);
    }
    else if (nrhs == 0) {
	    plhs[0] = mxCreateCellMatrix(3,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("Voltages"));
	    mxSetCell(plhs[0],2,mxCreateString("Frequency"));
	    if (nlhs > 1) { /* optional fields */
	        plhs[1] = mxCreateCellMatrix(4,1);
    	    mxSetCell(plhs[0],0,mxCreateString("Energy"));
            mxSetCell(plhs[1],1,mxCreateString("Phase"));
			mxSetCell(plhs[1],2,mxCreateString("SigPhi"));
			mxSetCell(plhs[1],3,mxCreateString("SigVV"));
	    }
    }
    else {
      mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
