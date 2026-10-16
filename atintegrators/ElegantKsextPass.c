/*
 * Canonical body integrator for elegant's KSEXT element.
 *
 * KSEXT is implemented by integrate_kick_multipole_ordn() in
 * elegant/src/multipole.c. Elegant exposes slopes at element boundaries,
 * but its body integrator uses canonical transverse momenta. PyAT already
 * uses these canonical coordinates, so the map can be applied directly.
 *
 * PyAT stores the sextupole coefficient as PolynomB[2] = K2/2, which
 * produces the same canonical kick as elegant's geometric K2 strength.
 * Fringe fields, radiation, and Elegant error-multipole files are excluded.
 *
 * Portions of the algorithm and formulas are based on Elegant source code:
 *   https://github.com/rtsoliday/elegant
 *
 * Elegant is copyright (c) 2002 University of Chicago. See
 * doc/licenseFiles/ElegantLicense.txt for the Elegant license terms.
 *
 * This file adapts the relevant Elegant method to AT/PyAT data structures,
 * canonical coordinates, aperture handling, and pass-method conventions.
 * Not all Elegant element features are implemented.
 */

#include "atelem.c"
#include "atlalib.c"

struct elem {
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    int IntegrationOrder;
    int ExpandHamiltonian;
    double Scaling;
    double *KickAngle;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *EApertures;
    double *RApertures;
};

static const double drift_frac_2[2] = {0.5, 0.5};
static const double kick_frac_2[2] = {1.0, 0.0};

static const double drift_frac_4[4] = {
    0.67560359597982889,
    -0.17560359597982883,
    -0.17560359597982883,
    0.67560359597982889
};
static const double kick_frac_4[4] = {
    1.3512071919596578,
    -1.7024143839193153,
    1.3512071919596578,
    0.0
};

static const double drift_frac_6[8] = {
    0.39225680523878,
    0.5100434119184585,
    -0.47105338540975655,
    0.0687531682525181,
    0.0687531682525181,
    -0.47105338540975655,
    0.5100434119184585,
    0.39225680523878
};
static const double kick_frac_6[8] = {
    0.784513610477560,
    0.235573213359357,
    -1.17767998417887,
    1.3151863206839063,
    -1.17767998417887,
    0.235573213359357,
    0.784513610477560,
    0.0
};

static int momenta_to_slopes(double px, double py, double delta,
                             double *xp, double *yp)
{
    double dp1 = 1.0 + delta;
    double pz2 = dp1 * dp1 - px * px - py * py;
    double pz;

    if (pz2 <= 0.0)
        return 0;

    pz = sqrt(pz2);
    *xp = px / pz;
    *yp = py / pz;
    return 1;
}

static void multipole_kick(double *r6, const struct elem *Elem,
                           double kick_length)
{
    double re_sum = Elem->PolynomB[Elem->MaxOrder];
    double im_sum = Elem->PolynomA[Elem->MaxOrder];
    int i;

    for (i = Elem->MaxOrder - 1; i >= 0; i--) {
        double re_next = re_sum * r6[x_] - im_sum * r6[y_]
                       + Elem->PolynomB[i];
        im_sum = im_sum * r6[x_] + re_sum * r6[y_]
               + Elem->PolynomA[i];
        re_sum = re_next;
    }

    r6[px_] -= kick_length * Elem->Scaling * re_sum;
    r6[py_] += kick_length * Elem->Scaling * im_sum;
}

static int body_map(double *r6, const struct elem *Elem)
{
    const double *drift_frac;
    const double *kick_frac;
    double slice_length = Elem->Length / Elem->NumIntSteps;
    double xp;
    double yp;
    double distance = 0.0;
    int substeps;
    int slice;
    int step;

    switch (Elem->IntegrationOrder) {
    case 2:
        drift_frac = drift_frac_2;
        kick_frac = kick_frac_2;
        substeps = 2;
        break;
    case 4:
        drift_frac = drift_frac_4;
        kick_frac = kick_frac_4;
        substeps = 4;
        break;
    default:
        drift_frac = drift_frac_6;
        kick_frac = kick_frac_6;
        substeps = 8;
        break;
    }

    if (!momenta_to_slopes(r6[px_], r6[py_], r6[delta_], &xp, &yp))
        return 0;

    for (slice = 0; slice < Elem->NumIntSteps; slice++) {
        for (step = 0; step < substeps; step++) {
            double drift_length = slice_length * drift_frac[step];

            r6[x_] += xp * drift_length;
            r6[y_] += yp * drift_length;
            if (Elem->ExpandHamiltonian)
                distance += drift_length
                          * (1.0 + 0.5 * (xp * xp + yp * yp));
            else
                distance += drift_length
                          * sqrt(1.0 + xp * xp + yp * yp);

            if (kick_frac[step] == 0.0)
                break;

            multipole_kick(r6, Elem,
                           slice_length * kick_frac[step]);
            if (Elem->KickAngle) {
                r6[px_] += Elem->KickAngle[0]
                         * kick_frac[step] / Elem->NumIntSteps;
                r6[py_] += Elem->KickAngle[1]
                         * kick_frac[step] / Elem->NumIntSteps;
            }

            if (!momenta_to_slopes(r6[px_], r6[py_], r6[delta_],
                                    &xp, &yp))
                return 0;
        }
    }

    r6[ct_] += distance - Elem->Length;
    return 1;
}

static void elegant_ksext(double *r_in, const struct elem *Elem,
                          int num_particles)
{
    int c;

    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) \
        default(none) shared(r_in, Elem, num_particles)
    for (c = 0; c < num_particles; c++) {
        double *r6 = r_in + 6 * c;

        if (!atIsNaN(r6[x_])) {
            if (Elem->T1)
                ATaddvv(r6, Elem->T1);
            if (Elem->R1)
                ATmultmv(r6, Elem->R1);
            if (Elem->RApertures)
                checkiflostRectangularAp(r6, Elem->RApertures);
            if (Elem->EApertures)
                checkiflostEllipticalAp(r6, Elem->EApertures);

            if (!atIsNaN(r6[x_]) && !body_map(r6, Elem))
                r6[x_] = atGetNaN();

            if (Elem->RApertures)
                checkiflostRectangularAp(r6, Elem->RApertures);
            if (Elem->EApertures)
                checkiflostEllipticalAp(r6, Elem->EApertures);
            if (Elem->R2)
                ATmultmv(r6, Elem->R2);
            if (Elem->T2)
                ATaddvv(r6, Elem->T2);
        }
    }
}

static struct elem *initialize_elem(const atElem *ElemData, struct elem *Elem)
{
    Elem->Length = atGetDouble(ElemData, "Length"); check_error();
    Elem->PolynomA = atGetDoubleArray(ElemData, "PolynomA"); check_error();
    Elem->PolynomB = atGetDoubleArray(ElemData, "PolynomB"); check_error();
    Elem->MaxOrder = atGetLong(ElemData, "MaxOrder"); check_error();
    Elem->NumIntSteps = atGetLong(ElemData, "NumIntSteps"); check_error();
    Elem->IntegrationOrder =
        atGetOptionalLong(ElemData, "IntegrationOrder", 4); check_error();
    Elem->ExpandHamiltonian =
        atGetOptionalLong(ElemData, "ExpandHamiltonian", 0); check_error();
    Elem->Scaling =
        atGetOptionalDouble(ElemData, "FieldScaling", 1.0); check_error();
    Elem->KickAngle =
        atGetOptionalDoubleArray(ElemData, "KickAngle"); check_error();
    Elem->R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
    Elem->R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
    Elem->T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
    Elem->T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();
    Elem->EApertures =
        atGetOptionalDoubleArray(ElemData, "EApertures"); check_error();
    Elem->RApertures =
        atGetOptionalDoubleArray(ElemData, "RApertures"); check_error();

    if (Elem->NumIntSteps < 1)
        atError("ElegantKsextPass: NumIntSteps must be positive");
    if (Elem->MaxOrder < 0)
        atError("ElegantKsextPass: MaxOrder must be nonnegative");
    if (Elem->IntegrationOrder != 2 && Elem->IntegrationOrder != 4
            && Elem->IntegrationOrder != 6)
        atError("ElegantKsextPass: IntegrationOrder must be 2, 4, or 6");
    check_error();
    return Elem;
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData, struct elem *Elem,
                                      double *r_in, int num_particles,
                                      struct parameters *Param)
{
    (void)Param;

    if (!Elem) {
        Elem = (struct elem *)atMalloc(sizeof(struct elem));
        Elem = initialize_elem(ElemData, Elem); check_error();
    }

    elegant_ksext(r_in, Elem, num_particles);
    return Elem;
}

MODULE_DEF(ElegantKsextPass)
#endif

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        struct elem Elem;
        double *r_in;
        int num_particles = mxGetN(prhs[1]);

        if (mxGetM(prhs[1]) != 6)
            mexErrMsgIdAndTxt("AT:WrongArg",
                              "Second argument must be a 6 x N matrix");
        initialize_elem(prhs[0], &Elem);
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        elegant_ksext(r_in, &Elem, num_particles);
    } else if (nrhs == 0) {
        plhs[0] = mxCreateCellMatrix(5, 1);
        mxSetCell(plhs[0], 0, mxCreateString("Length"));
        mxSetCell(plhs[0], 1, mxCreateString("PolynomA"));
        mxSetCell(plhs[0], 2, mxCreateString("PolynomB"));
        mxSetCell(plhs[0], 3, mxCreateString("MaxOrder"));
        mxSetCell(plhs[0], 4, mxCreateString("NumIntSteps"));
        if (nlhs > 1) {
            plhs[1] = mxCreateCellMatrix(10, 1);
            mxSetCell(plhs[1], 0, mxCreateString("IntegrationOrder"));
            mxSetCell(plhs[1], 1, mxCreateString("ExpandHamiltonian"));
            mxSetCell(plhs[1], 2, mxCreateString("FieldScaling"));
            mxSetCell(plhs[1], 3, mxCreateString("KickAngle"));
            mxSetCell(plhs[1], 4, mxCreateString("R1"));
            mxSetCell(plhs[1], 5, mxCreateString("R2"));
            mxSetCell(plhs[1], 6, mxCreateString("T1"));
            mxSetCell(plhs[1], 7, mxCreateString("T2"));
            mxSetCell(plhs[1], 8, mxCreateString("EApertures"));
            mxSetCell(plhs[1], 9, mxCreateString("RApertures"));
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif
