/*
 * Exact static corrector map corresponding to elegant's EKICKER element.
 *
 * The body map follows trackThroughExactCorrector() in
 * elegant/src/exactCorrector.c. Elegant evaluates the map in trace-space
 * slopes. PyAT coordinates are converted to slopes at the body entrance and
 * back to canonical momenta at the exit.
 *
 * Synchrotron radiation and steering-multipole error files are excluded.
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
    double Scaling;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *KickAngle;
};

static int momenta_to_slopes(const double *r6, double *xp, double *yp)
{
    double dp1 = 1.0 + r6[delta_];
    double pz2 = dp1 * dp1
               - r6[px_] * r6[px_]
               - r6[py_] * r6[py_];
    double pz;

    if (pz2 <= 0.0)
        return 0;

    pz = sqrt(pz2);
    *xp = r6[px_] / pz;
    *yp = r6[py_] / pz;
    return 1;
}

static void slopes_to_momenta(double *r6, double xp, double yp)
{
    double scale = (1.0 + r6[delta_])
                 / sqrt(1.0 + xp * xp + yp * yp);

    r6[px_] = xp * scale;
    r6[py_] = yp * scale;
}

static void rotate_coordinates(double *x, double *xp, double *y, double *yp,
                               double angle)
{
    double cosine = cos(angle);
    double sine = sin(angle);
    double x0 = *x;
    double xp0 = *xp;
    double y0 = *y;
    double yp0 = *yp;

    *x = x0 * cosine + y0 * sine;
    *y = -x0 * sine + y0 * cosine;
    *xp = xp0 * cosine + yp0 * sine;
    *yp = -xp0 * sine + yp0 * cosine;
}

static int exact_drift(double *r6, double length)
{
    double dp1 = 1.0 + r6[delta_];
    double pz2 = dp1 * dp1
               - r6[px_] * r6[px_]
               - r6[py_] * r6[py_];
    double normalized_length;

    if (pz2 <= 0.0)
        return 0;

    normalized_length = length / sqrt(pz2);
    r6[x_] += r6[px_] * normalized_length;
    r6[y_] += r6[py_] * normalized_length;
    r6[ct_] += dp1 * normalized_length - length;
    return 1;
}

static int exact_corrector(double *r6, double length,
                           double xkick, double ykick)
{
    double tan_xkick = tan(xkick);
    double tan_ykick = tan(ykick);
    double theta0 = atan(sqrt(tan_xkick * tan_xkick
                            + tan_ykick * tan_ykick));
    double tilt;
    double rho0;
    double x;
    double xp;
    double y;
    double yp;

    if (theta0 == 0.0)
        return exact_drift(r6, length);

    if (!momenta_to_slopes(r6, &xp, &yp))
        return 0;

    x = r6[x_];
    y = r6[y_];
    tilt = atan2(tan_ykick, tan_xkick);
    rho0 = fabs(length) / sin(theta0);
    if (length < 0.0)
        tilt += atan2(0.0, -1.0);

    rotate_coordinates(&x, &xp, &y, &yp, tilt);

    if (rho0 == 0.0 || length == 0.0) {
        xp = tan(theta0 / (1.0 + r6[delta_]) + atan(xp));
    } else {
        double alpha = atan(xp);
        double rho = rho0 * (1.0 + r6[delta_]);
        double arg = length / rho + sin(alpha);
        double theta;
        double dyf_denom;
        double dyf;
        double path_length;

        if (rho <= 0.0 || arg < -1.0 || arg > 1.0)
            return 0;

        theta = asin(arg) - alpha;
        dyf_denom = 1.0 + xp * xp - yp * yp;
        if (dyf_denom <= 0.0)
            return 0;
        dyf = yp / sqrt(dyf_denom);

        y += dyf * theta * rho;
        path_length = theta * rho * sqrt(1.0 + dyf * dyf);
        if (path_length < length)
            path_length = length;
        r6[ct_] += path_length - length;
        x += rho * (cos(alpha) - cos(theta + alpha));
        xp = tan(theta + alpha);
    }

    rotate_coordinates(&x, &xp, &y, &yp, -tilt);
    r6[x_] = x;
    r6[y_] = y;
    slopes_to_momenta(r6, xp, yp);
    return 1;
}

static void elegant_kick(double *r_in, const struct elem *Elem,
                         int num_particles)
{
    int c;

    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) \
        default(none) shared(r_in, Elem, num_particles)
    for (c = 0; c < num_particles; c++) {
        double *r6 = r_in + 6 * c;

        if (!atIsNaN(r6[x_])) {
            if (Elem->Scaling != 1.0)
                ATChangePRef(r6, Elem->Scaling);
            if (Elem->T1)
                ATaddvv(r6, Elem->T1);
            if (Elem->R1)
                ATmultmv(r6, Elem->R1);
            if (Elem->RApertures)
                checkiflostRectangularAp(r6, Elem->RApertures);
            if (Elem->EApertures)
                checkiflostEllipticalAp(r6, Elem->EApertures);

            if (!atIsNaN(r6[x_])
                    && !exact_corrector(r6, Elem->Length,
                                        Elem->KickAngle[0],
                                        Elem->KickAngle[1]))
                r6[x_] = atGetNaN();

            if (Elem->RApertures)
                checkiflostRectangularAp(r6, Elem->RApertures);
            if (Elem->EApertures)
                checkiflostEllipticalAp(r6, Elem->EApertures);
            if (Elem->R2)
                ATmultmv(r6, Elem->R2);
            if (Elem->T2)
                ATaddvv(r6, Elem->T2);
            if (Elem->Scaling != 1.0)
                ATChangePRef(r6, 1.0 / Elem->Scaling);
        }
    }
}

static struct elem *initialize_elem(const atElem *ElemData, struct elem *Elem)
{
    Elem->Length = atGetDouble(ElemData, "Length"); check_error();
    Elem->KickAngle = atGetDoubleArray(ElemData, "KickAngle"); check_error();
    Elem->Scaling =
        atGetOptionalDouble(ElemData, "FieldScaling", 1.0); check_error();
    Elem->R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
    Elem->R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
    Elem->T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
    Elem->T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();
    Elem->EApertures =
        atGetOptionalDoubleArray(ElemData, "EApertures"); check_error();
    Elem->RApertures =
        atGetOptionalDoubleArray(ElemData, "RApertures"); check_error();

    if (Elem->Scaling == 0.0)
        atError("ElegantEkickerPass: FieldScaling must be nonzero");
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

    elegant_kick(r_in, Elem, num_particles);
    return Elem;
}

MODULE_DEF(ElegantEkickerPass)
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
        elegant_kick(r_in, &Elem, num_particles);
    } else if (nrhs == 0) {
        plhs[0] = mxCreateCellMatrix(2, 1);
        mxSetCell(plhs[0], 0, mxCreateString("Length"));
        mxSetCell(plhs[0], 1, mxCreateString("KickAngle"));
        if (nlhs > 1) {
            plhs[1] = mxCreateCellMatrix(7, 1);
            mxSetCell(plhs[1], 0, mxCreateString("FieldScaling"));
            mxSetCell(plhs[1], 1, mxCreateString("R1"));
            mxSetCell(plhs[1], 2, mxCreateString("R2"));
            mxSetCell(plhs[1], 3, mxCreateString("T1"));
            mxSetCell(plhs[1], 4, mxCreateString("T2"));
            mxSetCell(plhs[1], 5, mxCreateString("EApertures"));
            mxSetCell(plhs[1], 6, mxCreateString("RApertures"));
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif
