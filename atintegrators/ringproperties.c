#include <math.h>
#include "atconstants.h"

static double atGetOptionalDoubleProp(const mxArray *obj, const char *fieldname, double default_value)
{
    mxArray *field=mxGetProperty(obj, 0, fieldname);
    return (field) ? mxGetScalar(field) : default_value;
}

static void atParticle(const mxArray *opts, double *energy, double *rest_energy, double *charge)
{
    const mxArray *part = mxGetField(opts, 0, "Particle");
    if (part) {
        if (mxIsClass(part, "atparticle")) {  /* OK */
            *rest_energy = atGetOptionalDoubleProp(part, "rest_energy", 0.0);
            *charge = atGetOptionalDoubleProp(part, "charge", -1.0);
            *energy = atGetOptionalDoubleProp(part, "energy", *energy);
        }
        else {                              /* particle is not a Particle object */
            mexErrMsgIdAndTxt("Atpass:WrongParameter","Particle must be an 'atparticle' object");
        }
    }
}

static void atProperties(const mxArray *opts, double *energy, double *rest_energy, double *charge)
{
    mxArray *field;
    if (!mxIsStruct(opts)) {
        mexErrMsgIdAndTxt("Atpass:WrongParameter","ring properties must be a struct");
    }
    atParticle(opts, energy, rest_energy, charge);
    field = mxGetField(opts, 0, "Energy");
    if (field) {
        double ener = mxGetScalar(field);
        if (ener > 0.0) {
            *energy = ener;
        }else{
            mexErrMsgIdAndTxt("Atpass:WrongParameter", "Energy must be positive. Check lattice or passmethod parameters.");
        }
    }
}

static void atParameters(
    const mxArray *opts,
    double *energy,
    double *rest_energy,
    double *charge,
    double *gamma0,
    double *beta0
)
{
    atProperties(opts, energy, rest_energy, charge);
    if (*rest_energy == 0.0) {
        *gamma0 = 1.0E-9 * *energy / __E0;
        *beta0 = 1.0;
    }
    else {
        double g0 = *energy / *rest_energy;
        double b0g0 = sqrt(g0*g0 - 1.0);
        *gamma0 = g0;
        *beta0 = b0g0/g0;
    }
}
