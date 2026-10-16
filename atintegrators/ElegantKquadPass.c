/*
 * Canonical body integrator for elegant's KQUAD element.
 *
 * KQUAD is implemented by integrate_kick_multipole_ordn() in
 * elegant/src/multipole.c. Elegant exposes slopes at element boundaries,
 * but its body integrator uses canonical transverse momenta. PyAT already
 * uses these canonical coordinates, so the map can be applied directly.
 *
 * Radiation and Elegant error-multipole files are excluded.
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
#include "quadfringe.c"

struct elem {
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    int IntegrationOrder;
    int ExpandHamiltonian;
    double Scaling;
    int FringeQuadEntrance;
    int FringeQuadExit;
    double *fringeIntM0;
    double *fringeIntP0;
    int edge1_order;
    int edge2_order;
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

struct fringe_matrix {
    double x00;
    double x01;
    double x10;
    double x11;
    double y00;
    double y01;
    double y10;
    double y11;
};

static void elegant_quad_partial_fringe_matrix(
    struct fringe_matrix *matrix, double k1, const double *fringe_int,
    int part)
{
    double j1x;
    double j2x;
    double j3x;
    double j1y;
    double j2y;
    double j3y;
    double k1_squared = k1 * k1;
    double exp_j1x;
    double exp_j1y;

    if (part == 1) {
        j1x = k1 * fringe_int[1]
            - 2.0 * k1_squared * fringe_int[3] / 3.0
            - k1_squared * fringe_int[0] * fringe_int[2] / 2.0;
        j2x = k1 * fringe_int[2];
        j3x = k1_squared
            * (fringe_int[2] + fringe_int[4]
               + fringe_int[0] * fringe_int[1]);

        j1y = -k1 * fringe_int[1]
            - 2.0 * k1_squared * fringe_int[3] / 3.0
            - k1_squared * fringe_int[0] * fringe_int[2] / 2.0;
        j2y = -j2x;
        j3y = j3x;
    } else {
        j1x = k1 * fringe_int[1]
            + k1_squared * fringe_int[0] * fringe_int[2] / 2.0;
        j2x = k1 * fringe_int[2];
        j3x = k1_squared
            * (fringe_int[4] - fringe_int[0] * fringe_int[1]);

        j1y = -k1 * fringe_int[1]
            + k1_squared * fringe_int[0] * fringe_int[2] / 2.0;
        j2y = -j2x;
        j3y = j3x;
    }

    exp_j1x = exp(j1x);
    matrix->x00 = exp_j1x;
    matrix->x01 = j2x / exp_j1x;
    matrix->x10 = exp_j1x * j3x;
    matrix->x11 = (1.0 + j2x * j3x) / exp_j1x;

    exp_j1y = exp(j1y);
    matrix->y00 = exp_j1y;
    matrix->y01 = j2y / exp_j1y;
    matrix->y10 = exp_j1y * j3y;
    matrix->y11 = (1.0 + j2y * j3y) / exp_j1y;
}

static void swap_double(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

static void orient_elegant_fringe_matrix(
    struct fringe_matrix *matrix, int reverse, int backtrack)
{
    if (reverse) {
        swap_double(&matrix->x00, &matrix->x11);
        swap_double(&matrix->y00, &matrix->y11);
    }
    if (backtrack) {
        swap_double(&matrix->x00, &matrix->x11);
        swap_double(&matrix->y00, &matrix->y11);
        matrix->x01 = -matrix->x01;
        matrix->x10 = -matrix->x10;
        matrix->y01 = -matrix->y01;
        matrix->y10 = -matrix->y10;
    }
}

static void apply_elegant_fringe_matrix(
    double *r6, const struct fringe_matrix *matrix)
{
    double x = matrix->x00 * r6[x_] + matrix->x01 * r6[px_];
    double px = matrix->x10 * r6[x_] + matrix->x11 * r6[px_];
    double y = matrix->y00 * r6[y_] + matrix->y01 * r6[py_];
    double py = matrix->y10 * r6[y_] + matrix->y11 * r6[py_];

    r6[x_] = x;
    r6[px_] = px;
    r6[y_] = y;
    r6[py_] = py;
}

static void elegant_quad_fringe(
    double *r6, double k1, const double *fringe_int_m,
    const double *fringe_int_p, int backtrack, int in_fringe,
    int higher_order)
{
    struct fringe_matrix matrix1;
    struct fringe_matrix matrix2;
    struct fringe_matrix tmp;
    double x;
    double px;
    double y;
    double py;
    double delta = r6[delta_];
    double a;
    double dx = 0.0;
    double dpx = 0.0;
    double dy = 0.0;
    double dpy = 0.0;
    double ds = 0.0;
    int reverse = in_fringe * (backtrack ? -1 : 1) == -1;
    int order_abs = higher_order < 0 ? -higher_order : higher_order;

    elegant_quad_partial_fringe_matrix(
        &matrix1, k1, fringe_int_m, 1);
    elegant_quad_partial_fringe_matrix(
        &matrix2, k1, fringe_int_p, 2);

    orient_elegant_fringe_matrix(&matrix1, reverse, backtrack);
    orient_elegant_fringe_matrix(&matrix2, reverse, backtrack);
    if (reverse) {
        tmp = matrix1;
        matrix1 = matrix2;
        matrix2 = tmp;
    }

    apply_elegant_fringe_matrix(r6, &matrix1);
    x = r6[x_];
    px = r6[px_];
    y = r6[y_];
    py = r6[py_];
    a = -in_fringe * k1 / (12.0 * (1.0 + delta));

    if (order_abs > 1) {
        double xpow[9];
        double ypow[9];
        double apow[4];
        int i;

        xpow[0] = 1.0;
        ypow[0] = 1.0;
        apow[0] = 1.0;

        if (order_abs > 2) {
            for (i = 1; i < 9; i++) {
                xpow[i] = xpow[i - 1] * x;
                ypow[i] = ypow[i - 1] * y;
            }
            for (i = 1; i < 4; i++)
                apow[i] = apow[i - 1] * a;

            dpx = (a * (
                12.0 * a
                    * (-4.0 * py * x * y
                       + px * (xpow[2] - 5.0 * ypow[2]))
                    * (xpow[2] - ypow[2])
                - 24.0
                    * (-2.0 * py * x * y
                       + px * (xpow[2] + ypow[2]))
                + 4.0 * apow[2]
                    * (-2.0 * py * x * y
                           * (3.0 * xpow[4]
                              + 50.0 * xpow[2] * ypow[2]
                              - 21.0 * ypow[4])
                       + px
                           * (xpow[6]
                              + 75.0 * xpow[4] * ypow[2]
                              - 105.0 * xpow[2] * ypow[4]
                              - 35.0 * ypow[6]))
                + 3.0 * apow[3]
                    * (-8.0 * py * x * y
                           * (xpow[6]
                              - 27.0 * xpow[4] * ypow[2]
                              + 83.0 * xpow[2] * ypow[4]
                              + 7.0 * ypow[6])
                       + px
                           * (xpow[8]
                              - 108.0 * xpow[6] * ypow[2]
                              + 830.0 * xpow[4] * ypow[4]
                              + 196.0 * xpow[2] * ypow[6]
                              + 105.0 * ypow[8])))) / 8.0;

            dpy = (a * (
                -8.0 * px * x * y
                    * (6.0
                       - 6.0 * a * (xpow[2] - ypow[2])
                       + apow[2]
                           * (21.0 * xpow[4]
                              - 50.0 * xpow[2] * ypow[2]
                              - 3.0 * ypow[4])
                       + 3.0 * apow[3]
                           * (7.0 * xpow[6]
                              + 83.0 * xpow[4] * ypow[2]
                              - 27.0 * xpow[2] * ypow[4]
                              + ypow[6]))
                + py
                    * (315.0 * apow[3] * xpow[8]
                       + 28.0 * apow[2] * xpow[6]
                           * (5.0 + 21.0 * a * ypow[2])
                       + 30.0 * a * xpow[4]
                           * (2.0
                              + 14.0 * a * ypow[2]
                              + 83.0 * apow[2] * ypow[4])
                       + ypow[2]
                           * (24.0
                              + 12.0 * a * ypow[2]
                              - 4.0 * apow[2] * ypow[4]
                              + 3.0 * apow[3] * ypow[6])
                       - 12.0 * xpow[2]
                           * (-2.0
                              + 6.0 * a * ypow[2]
                              + 25.0 * apow[2] * ypow[4]
                              + 27.0 * apow[3] * ypow[6])))) / 8.0;

            if (higher_order > 0) {
                double difference = xpow[2] - ypow[2];

                dx = (a * x * (
                    8.0 * (xpow[2] + 3.0 * ypow[2])
                    + 4.0 * apow[2]
                        * (5.0 * xpow[6]
                           + 21.0 * xpow[4] * ypow[2]
                           - 25.0 * xpow[2] * ypow[4]
                           - ypow[6])
                    + apow[3]
                        * (35.0 * xpow[8]
                           + 84.0 * xpow[6] * ypow[2]
                           + 498.0 * xpow[4] * ypow[4]
                           - 108.0 * xpow[2] * ypow[6]
                           + 3.0 * ypow[8])
                    + 12.0 * a * difference * difference)) / 8.0;
                dy = (a * y * (
                    -8.0 * (3.0 * xpow[2] + ypow[2])
                    + 4.0 * apow[2]
                        * (xpow[6]
                           + 25.0 * xpow[4] * ypow[2]
                           - 21.0 * xpow[2] * ypow[4]
                           - 5.0 * ypow[6])
                    + apow[3]
                        * (3.0 * xpow[8]
                           - 108.0 * xpow[6] * ypow[2]
                           + 498.0 * xpow[4] * ypow[4]
                           + 84.0 * xpow[2] * ypow[6]
                           + 35.0 * ypow[8])
                    + 12.0 * a * difference * difference)) / 8.0;
            }
        } else {
            for (i = 1; i < 4; i++) {
                xpow[i] = xpow[i - 1] * x;
                ypow[i] = ypow[i - 1] * y;
            }
            dpx = 3.0 * a
                * (-px * xpow[2] + 2.0 * py * x * y
                   - px * ypow[2]);
            dpy = 3.0 * a
                * (py * xpow[2] - 2.0 * px * x * y
                   + py * ypow[2]);
            if (higher_order > 0) {
                dx = a * (xpow[3] + 3.0 * x * ypow[2]);
                dy = a * (-3.0 * xpow[2] * y - ypow[3]);
            }
        }
        ds = (a / (1.0 + delta))
            * (3.0 * py * y * xpow[2]
               - px * xpow[3]
               - 3.0 * px * x * ypow[2]
               + py * ypow[3]);
    }

    r6[x_] = x + dx;
    r6[px_] = px + dpx;
    r6[y_] = y + dy;
    r6[py_] = py + dpy;
    apply_elegant_fringe_matrix(r6, &matrix2);
    r6[ct_] -= ds;
}

static void apply_quad_fringe(
    double *r6, const struct elem *Elem, int entrance)
{
    int method = entrance
        ? Elem->FringeQuadEntrance : Elem->FringeQuadExit;
    double k1;

    if (method == 0)
        return;

    /* Elegant applies FSE to the body, but not to quadFringe(). */
    k1 = Elem->PolynomB[1];
    if (k1 == 0.0)
        return;

    if (method == 1) {
        if (entrance)
            QuadFringePassP(r6, k1);
        else
            QuadFringePassN(r6, k1);
    } else if (method == 2) {
        if (entrance)
            linearQuadFringeElegantEntrance(
                r6, k1, Elem->fringeIntM0, Elem->fringeIntP0);
        else
            linearQuadFringeElegantExit(
                r6, k1, Elem->fringeIntM0, Elem->fringeIntP0);
    } else {
        elegant_quad_fringe(
            r6, k1, Elem->fringeIntM0, Elem->fringeIntP0,
            Elem->Length < 0.0, entrance ? -1 : 1,
            entrance ? Elem->edge1_order : Elem->edge2_order);
    }
}

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

static void elegant_kquad(double *r_in, const struct elem *Elem,
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

            if (!atIsNaN(r6[x_]))
                apply_quad_fringe(r6, Elem, 1);

            if (!atIsNaN(r6[x_]) && !body_map(r6, Elem))
                r6[x_] = atGetNaN();

            if (!atIsNaN(r6[x_]))
                apply_quad_fringe(r6, Elem, 0);

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
    Elem->FringeQuadEntrance =
        atGetOptionalLong(ElemData, "FringeQuadEntrance", 0);
    check_error();
    Elem->FringeQuadExit =
        atGetOptionalLong(ElemData, "FringeQuadExit", 0);
    check_error();
    Elem->fringeIntM0 =
        atGetOptionalDoubleArray(ElemData, "fringeIntM0");
    check_error();
    Elem->fringeIntP0 =
        atGetOptionalDoubleArray(ElemData, "fringeIntP0");
    check_error();
    Elem->edge1_order =
        atGetOptionalLong(ElemData, "edge1_order", 3); check_error();
    Elem->edge2_order =
        atGetOptionalLong(ElemData, "edge2_order", 3); check_error();
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
        atError("ElegantKquadPass: NumIntSteps must be positive");
    if (Elem->MaxOrder < 0)
        atError("ElegantKquadPass: MaxOrder must be nonnegative");
    if (Elem->IntegrationOrder != 2 && Elem->IntegrationOrder != 4
            && Elem->IntegrationOrder != 6)
        atError("ElegantKquadPass: IntegrationOrder must be 2, 4, or 6");
    if (Elem->FringeQuadEntrance < 0 || Elem->FringeQuadEntrance > 3
            || Elem->FringeQuadExit < 0 || Elem->FringeQuadExit > 3)
        atError("ElegantKquadPass: quadrupole fringe method must be 0 to 3");
    if ((Elem->FringeQuadEntrance || Elem->FringeQuadExit)
            && Elem->MaxOrder < 1)
        atError("ElegantKquadPass: quadrupole fringe requires MaxOrder >= 1");
    if ((Elem->FringeQuadEntrance >= 2 || Elem->FringeQuadExit >= 2)
            && (!Elem->fringeIntM0 || !Elem->fringeIntP0))
        atError("ElegantKquadPass: fringe methods 2 and 3 require "
                "fringeIntM0 and fringeIntP0");
    if (Elem->FringeQuadEntrance == 3
            && (Elem->edge1_order == 0 || Elem->edge1_order < -3
                || Elem->edge1_order > 3))
        atError("ElegantKquadPass: edge1_order must be +/-1, +/-2, or +/-3");
    if (Elem->FringeQuadExit == 3
            && (Elem->edge2_order == 0 || Elem->edge2_order < -3
                || Elem->edge2_order > 3))
        atError("ElegantKquadPass: edge2_order must be +/-1, +/-2, or +/-3");
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

    elegant_kquad(r_in, Elem, num_particles);
    return Elem;
}

MODULE_DEF(ElegantKquadPass)
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
        elegant_kquad(r_in, &Elem, num_particles);
    } else if (nrhs == 0) {
        plhs[0] = mxCreateCellMatrix(5, 1);
        mxSetCell(plhs[0], 0, mxCreateString("Length"));
        mxSetCell(plhs[0], 1, mxCreateString("PolynomA"));
        mxSetCell(plhs[0], 2, mxCreateString("PolynomB"));
        mxSetCell(plhs[0], 3, mxCreateString("MaxOrder"));
        mxSetCell(plhs[0], 4, mxCreateString("NumIntSteps"));
        if (nlhs > 1) {
            plhs[1] = mxCreateCellMatrix(16, 1);
            mxSetCell(plhs[1], 0, mxCreateString("IntegrationOrder"));
            mxSetCell(plhs[1], 1, mxCreateString("ExpandHamiltonian"));
            mxSetCell(plhs[1], 2, mxCreateString("FieldScaling"));
            mxSetCell(plhs[1], 3, mxCreateString("FringeQuadEntrance"));
            mxSetCell(plhs[1], 4, mxCreateString("FringeQuadExit"));
            mxSetCell(plhs[1], 5, mxCreateString("fringeIntM0"));
            mxSetCell(plhs[1], 6, mxCreateString("fringeIntP0"));
            mxSetCell(plhs[1], 7, mxCreateString("edge1_order"));
            mxSetCell(plhs[1], 8, mxCreateString("edge2_order"));
            mxSetCell(plhs[1], 9, mxCreateString("KickAngle"));
            mxSetCell(plhs[1], 10, mxCreateString("R1"));
            mxSetCell(plhs[1], 11, mxCreateString("R2"));
            mxSetCell(plhs[1], 12, mxCreateString("T1"));
            mxSetCell(plhs[1], 13, mxCreateString("T2"));
            mxSetCell(plhs[1], 14, mxCreateString("EApertures"));
            mxSetCell(plhs[1], 15, mxCreateString("RApertures"));
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif
