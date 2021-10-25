/* 
 *  RFCavityPass.c
 *  Accelerator Toolbox 
 *  22/09/2015
 *  Nicola Carmignani
 */

#include <math.h>
#include "atelem.c"
#include "atlalib.c"
#include "driftkickrad.c"

#define TWOPI  6.28318530717959
#define C0  	2.99792458e8 

struct elem 
{
  double Length;
  double NormV;
  double H_F;
  double Frequency;
  double TimeLag;
  double betgam0;
  double beta0;
  void (*RFKick)(double *r6, const struct elem* el, double tau_ref);
};

void non_relativistic_kick(double *r6, const struct elem *el, double tau_ref)
{
    double betgam = el.betgam0*(1.0 + r6[4]);
    double betai = betgam/sqrt(1.0 + betgami*betgami);
    double tau_rel = (r6[5]-el->TimeLag)/betai/C0;
    r6[4] += -el->NormV * sin(TWOPI*el->Frequency*(tau_rel + tau_ref));
}

void relativistic_kick(double *r6, const struct elem *el, double tau_ref)
{
    double tau_rel = (r6[5]-el->TimeLag)/C0;
    r6[4] += -el->NormV * sin(TWOPI*el->Frequency*(tau_rel + tau_ref));
}

void thin_cavity(double *r_in, const struct elem *el, double tau_ref, int num_particles)
{
    int c;
    for (c = 0; c < num_particles; c++) {
        double *r6 = r_in+c*6;
        if (!atIsNaN(r6[0]))
            /* Longitudinal momentum kick */
            el->RFKick(r6, el, tau_ref);
    }
}

void thick_cavity(double *r_in, const struct elem *el, double tau_ref, int num_particles)
{
    double halflength = el->Length/2;
    int c;
    for (c = 0; c < num_particles; c++) {
        double *r6 = r_in+c*6;
        if (!atIsNaN(r6[0]))  {
            /* Propagate through a drift equal to half cavity length */
            drift6(r6, halflength);
            /* Longitudinal momentum kick */
            el->RFKick(r6, el, tau_ref);
            /* Propagate through a drift equal to half cavity length */
            drift6(r6, halflength);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{
    int nturn=Param->nturn;
    double T0=Param->T0;
    if (!Elem) {
        double Length, Voltage, Energy, Frequency, TimeLag, gamma0;
        Length=atGetDouble(ElemData,"Length"); check_error();
        Voltage=atGetDouble(ElemData,"Voltage"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        Frequency=atGetDouble(ElemData,"Frequency"); check_error();
        TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0); check_error();
        gamma0 = Energy/param->rest_energy;
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->NormV=Voltage/Energy;
        Elem->Frequency=Frequency;
        Elem->H_F=round(Frequency*T0)/Frequency;
        Elem->TimeLag=TimeLag;
        if (isfinite(gamma0)) {
            Elem->betgam0 = sqrt(gamma0*gamma0 - 1.0);
            Elem->beta0 = Elem->betgam0/gamma0;
            Elem->RFKick = &nonrelativistic_kick;
        }
        else {
            Elem->RFKick = &relativistic_kick;
        }
    }
    double tau_ref = (T0-Elem->H_F)*nturn;
    if (Elem->Length == 0.0)
        thin_cavity(r_in, Elem, tau_ref, num_particles);
    else
        thick_cavity(r_in, Elem, tau_ref, num_particles);
    return Elem;
}

MODULE_DEF(RFCavityPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double energy=atGetDouble(ElemData,"Energy");
        double frequency=atGetDouble(ElemData,"Frequency");
        struct elem Elem;
        Elem.Length=atGetDouble(ElemData,"Length");
        Elem.normV=atGetDouble(ElemData,"Voltage")/energy;
        Elem.Frequency=frequency;
        Elem.TimeLag=atGetOptionalDouble(ElemData,"TimeLag",0);
        Elem.HarmNumber=round(frequency*T0)/frequency
        Elem.RFKick = &relativistic_kick;
        double T0=1.0/Frequency;      /* Does not matter since nturn == 0 */
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        double tau_ref = 0.0;
        if (Elem.Length == 0.0)
            thin_cavity(r_in, &Elem, tau_ref, num_particles);
        else
            thick_cavity(r_in, &Elem, tau_ref, num_particles);
    }
    else if (nrhs == 0) {
        /* return list of required fields */
        plhs[0] = mxCreateCellMatrix(4,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("Voltage"));
        mxSetCell(plhs[0],2,mxCreateString("Energy"));
        mxSetCell(plhs[0],3,mxCreateString("Frequency"));
        if (nlhs > 1) { /* optional fields */
            plhs[1] = mxCreateCellMatrix(1,1);
            mxSetCell(plhs[1],0,mxCreateString("TimeLag"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }

}
#endif /* MATLAB_MEX_FILE */
