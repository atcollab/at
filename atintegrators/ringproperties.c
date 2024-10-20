
static double atGetOptionalDoubleProp(const mxArray *obj, const char *fieldname, double default_value)
{
    mxArray *field=mxGetProperty(obj, 0, fieldname);
    return (field) ? mxGetScalar(field) : default_value;
}

static void atParticle(const mxArray *opts, double *rest_energy, double *charge)
{
    const mxArray *part = mxGetField(opts, 0, "Particle");
    if (part) {
        if (mxIsClass(part, "atparticle")) {  /* OK */
            *rest_energy = atGetOptionalDoubleProp(part, "rest_energy", 0.0);
            *charge = atGetOptionalDoubleProp(part, "charge", -1.0);
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
    field = mxGetField(opts, 0, "Energy");
    if (field) {
        double ener = mxGetScalar(field);
        if (ener != 0.0) *energy = ener;
        atParticle(opts, rest_energy, charge);
    }
}
