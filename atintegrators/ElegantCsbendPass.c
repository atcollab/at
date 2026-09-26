/*
 * Canonical body integrator for elegant's CSBEND element.
 *
 * This implements the non-expanded, non-radiative body map from
 * elegant/src/csbend.c. PyAT already uses the canonical transverse
 * coordinates required by the map, so no slope conversion is needed.
 *
 * The curvilinear field expansion supports normal and skew multipoles
 * through order 8 and total expansion order 10, matching CSBEND.
 * Edge effects are included for elegant's EDGE[12]_EFFECTS=4 model.
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

#define ELEGANT_MAX_EXPANSION_ORDER 10
#define ELEGANT_MAX_MULTIPOLE_ORDER 8

struct elem {
    double Length;
    double BendingAngle;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    int IntegrationOrder;
    int Nonlinear;
    int ExpansionOrder;
    double Scaling;
    double EntranceAngle;
    double ExitAngle;
    int FringeBendEntrance;
    int FringeBendExit;
    double FringeInt1;
    double FringeInt2;
    double FullGap;
    double H1;
    double H2;
    double FSE;
    double FSEDipole;
    double FSEQuadrupole;
    int FSECorrection;
    double FSECorrectionValue;
    double FSECorrectionPathError;
    double Fx[ELEGANT_MAX_EXPANSION_ORDER + 1]
             [ELEGANT_MAX_EXPANSION_ORDER + 1];
    double Fy[ELEGANT_MAX_EXPANSION_ORDER + 1]
             [ELEGANT_MAX_EXPANSION_ORDER + 1];
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

static double factorial(int order)
{
    double value = 1.0;
    int i;

    for (i = 2; i <= order; i++)
        value *= i;
    return value;
}

static double sqr(double x)
{
    return x * x;
}

static double ipow3(double x)
{
    return x * x * x;
}

static double safe_fse_denominator(double denominator)
{
    if (fabs(denominator) < 1.0e-12)
        return denominator < 0.0 ? -1.0e-12 : 1.0e-12;
    return denominator;
}

static double effective_fse(const struct elem *Elem)
{
    return Elem->FSE + Elem->FSEDipole
         + (Elem->FSECorrection ? Elem->FSECorrectionValue : 0.0);
}

static double actual_bending_radius(const struct elem *Elem)
{
    double rho0 = Elem->Length / Elem->BendingAngle;
    double scale = Elem->Scaling
                 * safe_fse_denominator(1.0 + effective_fse(Elem));

    return rho0 / scale;
}

static void dipole_fringe_khwang_rlindberg(double *r6, double rho,
                                           double in_fringe, double K1,
                                           double edge, double gap,
                                           double fint, double rhe)
{
    const double pi = 3.141592653589793238462643383279502884;
    double tan_edge, sin_edge, sec_edge, cos_edge;
    double cos3_edge, sec2_edge, tan2_edge, tan3_edge;
    double x0, px0, y0, py0, dp0;
    double x1, px1, y1, py1;
    double x2, px2, y2, py2;
    double x3, px3, y3, py3;
    double x4, px4, y4, py4;
    double x5, px5, y5, py5;
    double k0, k2, k3, k4, k5, k6;
    double t1, t2, rho2;

    (void)K1; /* Present in elegant's signature, unused in model 4. */

    if (rho == 0.0)
        return;
    if (gap == 0.0) {
        tan_edge = tan(edge);
        r6[px_] += tan_edge / rho * r6[x_];
        r6[py_] -= tan_edge / rho * r6[y_];
        return;
    }

    k0 = sqr(pi) / 6.0;
    k2 = fint;
    k3 = 1.0 / 6.0;
    k4 = -sqr(pi) / 3.0;
    k5 = 0.0;
    k6 = -1.0;
    rho2 = sqr(rho);

    x0 = r6[x_];
    px0 = r6[px_];
    y0 = r6[y_];
    py0 = r6[py_];
    dp0 = r6[delta_];

    cos_edge = cos(edge);
    sec_edge = 1.0 / cos_edge;
    tan_edge = tan(edge);
    sin_edge = sin(edge);

    sec2_edge = sqr(sec_edge);
    cos3_edge = ipow3(cos_edge);
    tan2_edge = sqr(tan_edge);
    tan3_edge = ipow3(tan_edge);

    if (in_fringe == -1.0) {
        x1 = x0;
        px1 = px0 + tan_edge / rho * x0
            + tan_edge / (2.0 * rho2 * (1.0 + dp0)) * sqr(y0)
            - tan3_edge / (rho2 * (1.0 + dp0)) * sqr(x0)
            - gap * k5 * sin_edge * rhe / (rho * cos3_edge) * x0
            + k6 * sec2_edge * rhe / (2.0 * rho)
            * (sqr(y0) - sqr(x0));
        y1 = y0;
        py1 = py0 - tan_edge / rho * y0
            + tan_edge / (rho2 * (1.0 + dp0)) * x0 * y0
            + gap * k2 * (1.0 + sqr(sin_edge))
            / (rho2 * (1.0 + dp0) * cos3_edge) * y0
            + 2.0 * k3 * (sqr(cos_edge) - 2.0)
            / (3.0 * gap * rho2 * cos3_edge) * ipow3(y0)
            + gap * k5 * sin_edge * rhe / (rho * cos3_edge) * y0
            + k6 * ipow3(sec_edge) * rhe / rho * x0 * y0;

        t1 = 1.0 + tan2_edge / (2.0 * rho * (1.0 + dp0)) * x1;
        x2 = x1 / t1;
        px2 = px1 * sqr(t1);
        y2 = y1;
        py2 = py1;

        t1 = sec2_edge / (rho * (1.0 + dp0));
        x3 = x2 + t1 / 2.0 * sqr(y2);
        px3 = px2;
        y3 = y2;
        py3 = py2 - t1 * y2 * px2;

        t1 = tan2_edge / (rho * (1.0 + dp0));
        t2 = exp(t1 * x3);
        x4 = x3;
        px4 = px3 - t1 * y3 * py3;
        y4 = y3 * t2;
        py4 = py3 / t2;

        t1 = sqr(gap) * k0 * sec2_edge / (rho * (1.0 + dp0));
        x5 = x4 - t1;
        px5 = px4 - t1 * tan_edge / rho
            + sqr(gap) * k4 * sqr(sin_edge) * rhe
            / (2.0 * rho * cos3_edge);
        y5 = y4;
        py5 = py4;
    } else {
        x1 = x0;
        px1 = px0 + tan_edge / rho * x0
            + tan3_edge / (2.0 * rho2 * (1.0 + dp0)) * sqr(y0)
            + tan3_edge / (2.0 * rho2 * (1.0 + dp0)) * sqr(x0)
            - gap * k5 * sin_edge * rhe / (rho * cos3_edge) * x0
            + k6 * sec2_edge * rhe / (2.0 * rho)
            * (sqr(y0) - sqr(x0));
        y1 = y0;
        py1 = py0 - tan_edge / rho * y0
            + tan3_edge / (rho2 * (1.0 + dp0)) * x0 * y0
            + gap * k2 * (1.0 + sqr(sin_edge))
            / (rho2 * (1.0 + dp0) * cos3_edge) * y0
            + 2.0 * k3 * (sqr(cos_edge) - 2.0)
            / (3.0 * gap * rho2 * cos3_edge) * ipow3(y0)
            + gap * k5 * sin_edge * rhe / (rho * cos3_edge) * y0
            + k6 * ipow3(sec_edge) * rhe / rho * x0 * y0;

        t1 = 1.0 - tan2_edge / (2.0 * rho * (1.0 + dp0)) * x1;
        x2 = x1 / t1;
        px2 = px1 * sqr(t1);
        y2 = y1;
        py2 = py1;

        t1 = sec2_edge / (rho * (1.0 + dp0));
        x3 = x2 - t1 / 2.0 * sqr(y2);
        px3 = px2;
        y3 = y2;
        py3 = py2 + t1 * y2 * px2;

        t1 = -tan2_edge / (rho * (1.0 + dp0));
        t2 = exp(t1 * x3);
        x4 = x3;
        px4 = px3 - t1 * y3 * py3;
        y4 = y3 * t2;
        py4 = py3 / t2;

        x5 = x4 + sqr(gap) * k0 * sec2_edge / (rho * (1.0 + dp0));
        px5 = px4 + sqr(gap) * k4 * sqr(sin_edge) * rhe
            / (2.0 * rho * cos3_edge);
        y5 = y4;
        py5 = py4;
    }

    r6[x_] = x5;
    r6[px_] = px5;
    r6[y_] = y5;
    r6[py_] = py5;
}

static void apply_lindberg_fringe(double *r6, const struct elem *Elem,
                                  double in_fringe)
{
    double rho_actual = actual_bending_radius(Elem);
    double K1 = Elem->MaxOrder >= 1 ? Elem->PolynomB[1] : 0.0;

    if (in_fringe == -1.0) {
        if (Elem->FringeBendEntrance == 4)
            dipole_fringe_khwang_rlindberg(r6, rho_actual, -1.0, K1,
                                           Elem->EntranceAngle,
                                           Elem->FullGap, Elem->FringeInt1,
                                           Elem->H1);
    } else {
        if (Elem->FringeBendExit == 4)
            dipole_fringe_khwang_rlindberg(r6, rho_actual, 1.0, K1,
                                           Elem->ExitAngle,
                                           Elem->FullGap, Elem->FringeInt2,
                                           Elem->H2);
    }
}

static void compute_field_coefficients(struct elem *Elem)
{
    double b[ELEGANT_MAX_MULTIPOLE_ORDER + 1] = {0.0};
    double c[ELEGANT_MAX_MULTIPOLE_ORDER + 1] = {0.0};
    double *b_formula;
    double *c_formula;
    double h[20];
    double rho0 = Elem->Length / Elem->BendingAngle;
    int i;
    int j;

    for (i = 0; i <= ELEGANT_MAX_EXPANSION_ORDER; i++)
        for (j = 0; j <= ELEGANT_MAX_EXPANSION_ORDER; j++)
            Elem->Fx[i][j] = Elem->Fy[i][j] = 0.0;

    for (i = 0; i <= Elem->MaxOrder; i++) {
        double scale = factorial(i) * rho0;
        b[i] = Elem->PolynomB[i] * scale;
        c[i] = Elem->PolynomA[i] * scale;
    }

    {
        double denominator =
            safe_fse_denominator(1.0 + Elem->FSE + Elem->FSEDipole);

        b[1] *= (1.0 + Elem->FSE + Elem->FSEQuadrupole) / denominator;
        c[1] *= (1.0 + Elem->FSE + Elem->FSEQuadrupole) / denominator;
        b[2] *= (1.0 + Elem->FSE) / denominator;
        c[2] *= (1.0 + Elem->FSE) / denominator;
    }

    /*
     * elegant defines b[0] as a fractional reduction of the design dipole
     * field, whereas AT PolynomB[0] is an additive normal dipole field.
     */
    b[0] = -b[0];

    if (Elem->ExpansionOrder == 0) {
        for (i = ELEGANT_MAX_MULTIPOLE_ORDER; i >= 0; i--)
            if (b[i] != 0.0 || c[i] != 0.0)
                break;
        Elem->ExpansionOrder = i + 2;
        if (Elem->ExpansionOrder < 4)
            Elem->ExpansionOrder = 4;
    }

    h[0] = 1.0;
    for (i = 1; i < 20; i++)
        h[i] = h[i - 1] / rho0;

    Elem->Fx[0][0] = c[0];
    Elem->Fy[0][0] = 1.0 - b[0];

    b_formula = b + 1;
    c_formula = c + 1;
    Elem->Fx[0][1] = b_formula[0];
    Elem->Fy[1][0] = b_formula[0];
    Elem->Fy[0][1] = c_formula[0];
    Elem->Fx[1][0] = -c_formula[0];

    /*
     * Source-generated CSBEND expressions. b_formula and c_formula use
     * elegant's quadrupole-first indexing.
     */
    if (Elem->Nonlinear) {
        Elem->Fy[0][2] = -(h[1] * b_formula[0]) / 2 - b_formula[1] / 2;
        Elem->Fy[0][3] = (h[2] * c_formula[0]) / 6 - (h[1] * c_formula[1]) / 3 - c_formula[2] / 6;
        Elem->Fy[0][4] = (h[3] * b_formula[0]) / 24 - (h[2] * b_formula[1]) / 24 + (h[1] * b_formula[2]) / 12 + b_formula[3] / 24;
        Elem->Fy[0][5] = (-3 * h[4] * c_formula[0]) / 40 + (h[3] * c_formula[1]) / 20 - (h[2] * c_formula[2]) / 40 + (h[1] * c_formula[3]) / 40 + c_formula[4] / 120;
        Elem->Fy[0][6] = -(h[5] * b_formula[0]) / 80 + (h[4] * b_formula[1]) / 80 - (h[3] * b_formula[2]) / 120 + (h[2] * b_formula[3]) / 240 - (h[1] * b_formula[4]) / 240 - b_formula[5] / 720;
        Elem->Fy[0][7] = (5 * h[6] * c_formula[0]) / 112 - (h[5] * c_formula[1]) / 40 + (17 * h[4] * c_formula[2]) / 1680 - (h[3] * c_formula[3]) / 280 + (h[2] * c_formula[4]) / 840 -
                      (h[1] * c_formula[5]) / 1260 - c_formula[6] / 5040;
        Elem->Fy[0][8] = (5 * h[7] * b_formula[0]) / 896 - (5 * h[6] * b_formula[1]) / 896 + (h[5] * b_formula[2]) / 320 - (17 * h[4] * b_formula[3]) / 13440 + (h[3] * b_formula[4]) / 2240 -
                      (h[2] * b_formula[5]) / 6720 + (h[1] * b_formula[6]) / 10080 + b_formula[7] / 40320;
        Elem->Fy[0][9] = (-35 * h[8] * c_formula[0]) / 1152 + (65 * h[7] * c_formula[1]) / 4032 - (145 * h[6] * c_formula[2]) / 24192 + (43 * h[5] * c_formula[3]) / 24192 -
                      (11 * h[4] * c_formula[4]) / 24192 + (h[3] * c_formula[5]) / 9072 - (h[2] * c_formula[6]) / 36288 + (h[1] * c_formula[7]) / 72576;
        Elem->Fy[0][10] = (-7 * h[9] * b_formula[0]) / 2304 + (7 * h[8] * b_formula[1]) / 2304 - (13 * h[7] * b_formula[2]) / 8064 + (29 * h[6] * b_formula[3]) / 48384 -
                       (43 * h[5] * b_formula[4]) / 241920 + (11 * h[4] * b_formula[5]) / 241920 - (h[3] * b_formula[6]) / 90720 + (h[2] * b_formula[7]) / 362880;
        Elem->Fy[1][1] = h[1] * c_formula[0] + c_formula[1];
        Elem->Fy[1][2] = (h[2] * b_formula[0]) / 2 - (h[1] * b_formula[1]) / 2 - b_formula[2] / 2;
        Elem->Fy[1][3] = -(h[3] * c_formula[0]) / 2 + (h[2] * c_formula[1]) / 2 - (h[1] * c_formula[2]) / 3 - c_formula[3] / 6;
        Elem->Fy[1][4] = -(h[4] * b_formula[0]) / 8 + (h[3] * b_formula[1]) / 8 - (h[2] * b_formula[2]) / 8 + (h[1] * b_formula[3]) / 12 + b_formula[4] / 24;
        Elem->Fy[1][5] = (3 * h[5] * c_formula[0]) / 8 - (9 * h[4] * c_formula[1]) / 40 + (h[3] * c_formula[2]) / 10 - (h[2] * c_formula[3]) / 20 + (h[1] * c_formula[4]) / 40 + c_formula[5] / 120;
        Elem->Fy[1][6] = (h[6] * b_formula[0]) / 16 - (h[5] * b_formula[1]) / 16 + (3 * h[4] * b_formula[2]) / 80 - (h[3] * b_formula[3]) / 60 + (h[2] * b_formula[4]) / 120 - (h[1] * b_formula[5]) / 240 -
                      b_formula[6] / 720;
        Elem->Fy[1][7] = (-5 * h[7] * c_formula[0]) / 16 + (19 * h[6] * c_formula[1]) / 112 - (11 * h[5] * c_formula[2]) / 168 + (h[4] * c_formula[3]) / 48 - (h[3] * c_formula[4]) / 168 +
                      (h[2] * c_formula[5]) / 504 - (h[1] * c_formula[6]) / 1260 - c_formula[7] / 5040;
        Elem->Fy[1][8] = (-5 * h[8] * b_formula[0]) / 128 + (5 * h[7] * b_formula[1]) / 128 - (19 * h[6] * b_formula[2]) / 896 + (11 * h[5] * b_formula[3]) / 1344 - (h[4] * b_formula[4]) / 384 +
                      (h[3] * b_formula[5]) / 1344 - (h[2] * b_formula[6]) / 4032 + (h[1] * b_formula[7]) / 10080;
        Elem->Fy[1][9] = (35 * h[9] * c_formula[0]) / 128 - (55 * h[8] * c_formula[1]) / 384 + (5 * h[7] * c_formula[2]) / 96 - (5 * h[6] * c_formula[3]) / 336 + (29 * h[5] * c_formula[4]) / 8064 -
                      (19 * h[4] * c_formula[5]) / 24192 + (h[3] * c_formula[6]) / 6048 - (h[2] * c_formula[7]) / 24192;
        Elem->Fy[1][10] = (7 * h[10] * b_formula[0]) / 256 - (7 * h[9] * b_formula[1]) / 256 + (11 * h[8] * b_formula[2]) / 768 - (h[7] * b_formula[3]) / 192 + (h[6] * b_formula[4]) / 672 -
                       (29 * h[5] * b_formula[5]) / 80640 + (19 * h[4] * b_formula[6]) / 241920 - (h[3] * b_formula[7]) / 60480;
        Elem->Fy[2][0] = b_formula[1] / 2;
        Elem->Fy[2][1] = -(h[2] * c_formula[0]) + (h[1] * c_formula[1]) / 2 + c_formula[2] / 2;
        Elem->Fy[2][2] = -(h[3] * b_formula[0]) / 2 + (h[2] * b_formula[1]) / 2 - (h[1] * b_formula[2]) / 4 - b_formula[3] / 4;
        Elem->Fy[2][3] = h[4] * c_formula[0] - (3 * h[3] * c_formula[1]) / 4 + (5 * h[2] * c_formula[2]) / 12 - (h[1] * c_formula[3]) / 6 - c_formula[4] / 12;
        Elem->Fy[2][4] = (h[5] * b_formula[0]) / 4 - (h[4] * b_formula[1]) / 4 + (3 * h[3] * b_formula[2]) / 16 - (5 * h[2] * b_formula[3]) / 48 + (h[1] * b_formula[4]) / 24 + b_formula[5] / 48;
        Elem->Fy[2][5] = (-9 * h[6] * c_formula[0]) / 8 + (51 * h[5] * c_formula[1]) / 80 - (21 * h[4] * c_formula[2]) / 80 + (h[3] * c_formula[3]) / 10 - (3 * h[2] * c_formula[4]) / 80 +
                      (h[1] * c_formula[5]) / 80 + c_formula[6] / 240;
        Elem->Fy[2][6] = (-3 * h[7] * b_formula[0]) / 16 + (3 * h[6] * b_formula[1]) / 16 - (17 * h[5] * b_formula[2]) / 160 + (7 * h[4] * b_formula[3]) / 160 - (h[3] * b_formula[4]) / 60 +
                      (h[2] * b_formula[5]) / 160 - (h[1] * b_formula[6]) / 480 - b_formula[7] / 1440;
        Elem->Fy[2][7] = (5 * h[8] * c_formula[0]) / 4 - (149 * h[7] * c_formula[1]) / 224 + (167 * h[6] * c_formula[2]) / 672 - (25 * h[5] * c_formula[3]) / 336 + (13 * h[4] * c_formula[4]) / 672 -
                      (5 * h[3] * c_formula[5]) / 1008 + (h[2] * c_formula[6]) / 720 - (h[1] * c_formula[7]) / 2520;
        Elem->Fy[2][8] = (5 * h[9] * b_formula[0]) / 32 - (5 * h[8] * b_formula[1]) / 32 + (149 * h[7] * b_formula[2]) / 1792 - (167 * h[6] * b_formula[3]) / 5376 +
                      (25 * h[5] * b_formula[4]) / 2688 - (13 * h[4] * b_formula[5]) / 5376 + (5 * h[3] * b_formula[6]) / 8064 - (h[2] * b_formula[7]) / 5760;
        Elem->Fy[2][9] = (-175 * h[10] * c_formula[0]) / 128 + (545 * h[9] * c_formula[1]) / 768 - (65 * h[8] * c_formula[2]) / 256 + (95 * h[7] * c_formula[3]) / 1344 -
                      (265 * h[6] * c_formula[4]) / 16128 + (163 * h[5] * c_formula[5]) / 48384 - (31 * h[4] * c_formula[6]) / 48384 + (h[3] * c_formula[7]) / 8064;
        Elem->Fy[2][10] = (-35 * h[11] * b_formula[0]) / 256 + (35 * h[10] * b_formula[1]) / 256 - (109 * h[9] * b_formula[2]) / 1536 + (13 * h[8] * b_formula[3]) / 512 -
                       (19 * h[7] * b_formula[4]) / 2688 + (53 * h[6] * b_formula[5]) / 32256 - (163 * h[5] * b_formula[6]) / 483840 + (31 * h[4] * b_formula[7]) / 483840;
        Elem->Fy[3][0] = b_formula[2] / 6;
        Elem->Fy[3][1] = h[3] * c_formula[0] - (h[2] * c_formula[1]) / 2 + (h[1] * c_formula[2]) / 6 + c_formula[3] / 6;
        Elem->Fy[3][2] = (h[4] * b_formula[0]) / 2 - (h[3] * b_formula[1]) / 2 + (h[2] * b_formula[2]) / 4 - (h[1] * b_formula[3]) / 12 - b_formula[4] / 12;
        Elem->Fy[3][3] = (-5 * h[5] * c_formula[0]) / 3 + (13 * h[4] * c_formula[1]) / 12 - (19 * h[3] * c_formula[2]) / 36 + (7 * h[2] * c_formula[3]) / 36 - (h[1] * c_formula[4]) / 18 - c_formula[5] / 36;
        Elem->Fy[3][4] = (-5 * h[6] * b_formula[0]) / 12 + (5 * h[5] * b_formula[1]) / 12 - (13 * h[4] * b_formula[2]) / 48 + (19 * h[3] * b_formula[3]) / 144 - (7 * h[2] * b_formula[4]) / 144 +
                      (h[1] * b_formula[5]) / 72 + b_formula[6] / 144;
        Elem->Fy[3][5] = (21 * h[7] * c_formula[0]) / 8 - (23 * h[6] * c_formula[1]) / 16 + (9 * h[5] * c_formula[2]) / 16 - (3 * h[4] * c_formula[3]) / 16 + (7 * h[3] * c_formula[4]) / 120 -
                      (h[2] * c_formula[5]) / 60 + (h[1] * c_formula[6]) / 240 + c_formula[7] / 720;
        Elem->Fy[3][6] = (7 * h[8] * b_formula[0]) / 16 - (7 * h[7] * b_formula[1]) / 16 + (23 * h[6] * b_formula[2]) / 96 - (3 * h[5] * b_formula[3]) / 32 + (h[4] * b_formula[4]) / 32 -
                      (7 * h[3] * b_formula[5]) / 720 + (h[2] * b_formula[6]) / 360 - (h[1] * b_formula[7]) / 1440;
        Elem->Fy[3][7] = (-15 * h[9] * c_formula[0]) / 4 + (63 * h[8] * c_formula[1]) / 32 - (23 * h[7] * c_formula[2]) / 32 + (139 * h[6] * c_formula[3]) / 672 - (17 * h[5] * c_formula[4]) / 336 +
                      (23 * h[4] * c_formula[5]) / 2016 - (13 * h[3] * c_formula[6]) / 5040 + (h[2] * c_formula[7]) / 1680;
        Elem->Fy[3][8] = (-15 * h[10] * b_formula[0]) / 32 + (15 * h[9] * b_formula[1]) / 32 - (63 * h[8] * b_formula[2]) / 256 + (23 * h[7] * b_formula[3]) / 256 -
                      (139 * h[6] * b_formula[4]) / 5376 + (17 * h[5] * b_formula[5]) / 2688 - (23 * h[4] * b_formula[6]) / 16128 + (13 * h[3] * b_formula[7]) / 40320;
        Elem->Fy[3][9] = (1925 * h[11] * c_formula[0]) / 384 - (1985 * h[10] * c_formula[1]) / 768 + (2105 * h[9] * c_formula[2]) / 2304 - (575 * h[8] * c_formula[3]) / 2304 +
                      (65 * h[7] * c_formula[4]) / 1152 - (115 * h[6] * c_formula[5]) / 10368 + (41 * h[5] * c_formula[6]) / 20736 - (7 * h[4] * c_formula[7]) / 20736;
        Elem->Fy[3][10] = (385 * h[12] * b_formula[0]) / 768 - (385 * h[11] * b_formula[1]) / 768 + (397 * h[10] * b_formula[2]) / 1536 - (421 * h[9] * b_formula[3]) / 4608 +
                       (115 * h[8] * b_formula[4]) / 4608 - (13 * h[7] * b_formula[5]) / 2304 + (23 * h[6] * b_formula[6]) / 20736 - (41 * h[5] * b_formula[7]) / 207360;
        Elem->Fy[4][0] = b_formula[3] / 24;
        Elem->Fy[4][1] = -(h[4] * c_formula[0]) + (h[3] * c_formula[1]) / 2 - (h[2] * c_formula[2]) / 6 + (h[1] * c_formula[3]) / 24 + c_formula[4] / 24;
        Elem->Fy[4][2] = -(h[5] * b_formula[0]) / 2 + (h[4] * b_formula[1]) / 2 - (h[3] * b_formula[2]) / 4 + (h[2] * b_formula[3]) / 12 - (h[1] * b_formula[4]) / 48 - b_formula[5] / 48;
        Elem->Fy[4][3] = (5 * h[6] * c_formula[0]) / 2 - (3 * h[5] * c_formula[1]) / 2 + (2 * h[4] * c_formula[2]) / 3 - (11 * h[3] * c_formula[3]) / 48 + (h[2] * c_formula[4]) / 16 - (h[1] * c_formula[5]) / 72 -
                      c_formula[6] / 144;
        Elem->Fy[4][4] = (5 * h[7] * b_formula[0]) / 8 - (5 * h[6] * b_formula[1]) / 8 + (3 * h[5] * b_formula[2]) / 8 - (h[4] * b_formula[3]) / 6 + (11 * h[3] * b_formula[4]) / 192 -
                      (h[2] * b_formula[5]) / 64 + (h[1] * b_formula[6]) / 288 + b_formula[7] / 576;
        Elem->Fy[4][5] = (-21 * h[8] * c_formula[0]) / 4 + (45 * h[7] * c_formula[1]) / 16 - (17 * h[6] * c_formula[2]) / 16 + (21 * h[5] * c_formula[3]) / 64 - (29 * h[4] * c_formula[4]) / 320 +
                      (11 * h[3] * c_formula[5]) / 480 - (h[2] * c_formula[6]) / 192 + (h[1] * c_formula[7]) / 960;
        Elem->Fy[4][6] = (-7 * h[9] * b_formula[0]) / 8 + (7 * h[8] * b_formula[1]) / 8 - (15 * h[7] * b_formula[2]) / 32 + (17 * h[6] * b_formula[3]) / 96 - (7 * h[5] * b_formula[4]) / 128 +
                      (29 * h[4] * b_formula[5]) / 1920 - (11 * h[3] * b_formula[6]) / 2880 + (h[2] * b_formula[7]) / 1152;
        Elem->Fy[4][7] = (75 * h[10] * c_formula[0]) / 8 - (39 * h[9] * c_formula[1]) / 8 + (7 * h[8] * c_formula[2]) / 4 - (439 * h[7] * c_formula[3]) / 896 + (103 * h[6] * c_formula[4]) / 896 -
                      (97 * h[5] * c_formula[5]) / 4032 + (193 * h[4] * c_formula[6]) / 40320 - (19 * h[3] * c_formula[7]) / 20160;
        Elem->Fy[4][8] = (75 * h[11] * b_formula[0]) / 64 - (75 * h[10] * b_formula[1]) / 64 + (39 * h[9] * b_formula[2]) / 64 - (7 * h[8] * b_formula[3]) / 32 + (439 * h[7] * b_formula[4]) / 7168 -
                      (103 * h[6] * b_formula[5]) / 7168 + (97 * h[5] * b_formula[6]) / 32256 - (193 * h[4] * b_formula[7]) / 322560;
        Elem->Fy[4][9] = (-1925 * h[12] * c_formula[0]) / 128 + (1975 * h[11] * c_formula[1]) / 256 - (2075 * h[10] * c_formula[2]) / 768 + (745 * h[9] * c_formula[3]) / 1024 -
                      (165 * h[8] * c_formula[4]) / 1024 + (425 * h[7] * c_formula[5]) / 13824 - (145 * h[6] * c_formula[6]) / 27648 + (23 * h[5] * c_formula[7]) / 27648;
        Elem->Fy[4][10] = (-385 * h[13] * b_formula[0]) / 256 + (385 * h[12] * b_formula[1]) / 256 - (395 * h[11] * b_formula[2]) / 512 + (415 * h[10] * b_formula[3]) / 1536 -
                       (149 * h[9] * b_formula[4]) / 2048 + (33 * h[8] * b_formula[5]) / 2048 - (85 * h[7] * b_formula[6]) / 27648 + (29 * h[6] * b_formula[7]) / 55296;
        Elem->Fy[5][0] = b_formula[4] / 120;
        Elem->Fy[5][1] = h[5] * c_formula[0] - (h[4] * c_formula[1]) / 2 + (h[3] * c_formula[2]) / 6 - (h[2] * c_formula[3]) / 24 + (h[1] * c_formula[4]) / 120 + c_formula[5] / 120;
        Elem->Fy[5][2] = (h[6] * b_formula[0]) / 2 - (h[5] * b_formula[1]) / 2 + (h[4] * b_formula[2]) / 4 - (h[3] * b_formula[3]) / 12 + (h[2] * b_formula[4]) / 48 - (h[1] * b_formula[5]) / 240 -
                      b_formula[6] / 240;
        Elem->Fy[5][3] = (-7 * h[7] * c_formula[0]) / 2 + 2 * h[6] * c_formula[1] - (5 * h[5] * c_formula[2]) / 6 + (13 * h[4] * c_formula[3]) / 48 - (17 * h[3] * c_formula[4]) / 240 +
                      (11 * h[2] * c_formula[5]) / 720 - (h[1] * c_formula[6]) / 360 - c_formula[7] / 720;
        Elem->Fy[5][4] = (-7 * h[8] * b_formula[0]) / 8 + (7 * h[7] * b_formula[1]) / 8 - (h[6] * b_formula[2]) / 2 + (5 * h[5] * b_formula[3]) / 24 - (13 * h[4] * b_formula[4]) / 192 +
                      (17 * h[3] * b_formula[5]) / 960 - (11 * h[2] * b_formula[6]) / 2880 + (h[1] * b_formula[7]) / 1440;
        Elem->Fy[5][5] = (189 * h[9] * c_formula[0]) / 20 - (399 * h[8] * c_formula[1]) / 80 + (147 * h[7] * c_formula[2]) / 80 - (173 * h[6] * c_formula[3]) / 320 +
                      (221 * h[5] * c_formula[4]) / 1600 - (51 * h[4] * c_formula[5]) / 1600 + (h[3] * c_formula[6]) / 150 - (h[2] * c_formula[7]) / 800;
        Elem->Fy[5][6] = (63 * h[10] * b_formula[0]) / 40 - (63 * h[9] * b_formula[1]) / 40 + (133 * h[8] * b_formula[2]) / 160 - (49 * h[7] * b_formula[3]) / 160 +
                      (173 * h[6] * b_formula[4]) / 1920 - (221 * h[5] * b_formula[5]) / 9600 + (17 * h[4] * b_formula[6]) / 3200 - (h[3] * b_formula[7]) / 900;
        Elem->Fy[5][7] = (-165 * h[11] * c_formula[0]) / 8 + (213 * h[10] * c_formula[1]) / 20 - (151 * h[9] * c_formula[2]) / 40 + (663 * h[8] * c_formula[3]) / 640 -
                      (151 * h[7] * c_formula[4]) / 640 + (271 * h[6] * c_formula[5]) / 5760 - (871 * h[5] * c_formula[6]) / 100800 + (307 * h[4] * c_formula[7]) / 201600;
        Elem->Fy[5][8] = (-165 * h[12] * b_formula[0]) / 64 + (165 * h[11] * b_formula[1]) / 64 - (213 * h[10] * b_formula[2]) / 160 + (151 * h[9] * b_formula[3]) / 320 -
                      (663 * h[8] * b_formula[4]) / 5120 + (151 * h[7] * b_formula[5]) / 5120 - (271 * h[6] * b_formula[6]) / 46080 + (871 * h[5] * b_formula[7]) / 806400;
        Elem->Fy[5][9] = (-13299 * h[13] * c_formula[0]) / 128 + (13189 * h[12] * c_formula[1]) / 256 - (4323 * h[11] * c_formula[2]) / 256 + (4207 * h[10] * c_formula[3]) / 1024 -
                      (12109 * h[9] * c_formula[4]) / 15360 + (17051 * h[8] * c_formula[5]) / 138240 - (1927 * h[7] * c_formula[6]) / 120960 + (403 * h[6] * c_formula[7]) / 241920;
        Elem->Fy[5][10] = (-13299 * h[14] * b_formula[0]) / 1280 + (13299 * h[13] * b_formula[1]) / 1280 - (13189 * h[12] * b_formula[2]) / 2560 +
                       (4323 * h[11] * b_formula[3]) / 2560 - (4207 * h[10] * b_formula[4]) / 10240 + (12109 * h[9] * b_formula[5]) / 153600 - (17051 * h[8] * b_formula[6]) / 1382400 +
                       (1927 * h[7] * b_formula[7]) / 1209600;
        Elem->Fy[6][0] = b_formula[5] / 720;
        Elem->Fy[6][1] = -(h[6] * c_formula[0]) + (h[5] * c_formula[1]) / 2 - (h[4] * c_formula[2]) / 6 + (h[3] * c_formula[3]) / 24 - (h[2] * c_formula[4]) / 120 + (h[1] * c_formula[5]) / 720 +
                      c_formula[6] / 720;
        Elem->Fy[6][2] = -(h[7] * b_formula[0]) / 2 + (h[6] * b_formula[1]) / 2 - (h[5] * b_formula[2]) / 4 + (h[4] * b_formula[3]) / 12 - (h[3] * b_formula[4]) / 48 + (h[2] * b_formula[5]) / 240 -
                      (h[1] * b_formula[6]) / 1440 - b_formula[7] / 1440;
        Elem->Fy[6][3] = (14 * h[8] * c_formula[0]) / 3 - (31 * h[7] * c_formula[1]) / 12 + (37 * h[6] * c_formula[2]) / 36 - (23 * h[5] * c_formula[3]) / 72 + (29 * h[4] * c_formula[4]) / 360 -
                      (73 * h[3] * c_formula[5]) / 4320 + (13 * h[2] * c_formula[6]) / 4320 - (h[1] * c_formula[7]) / 2160;
        Elem->Fy[6][4] = (7 * h[9] * b_formula[0]) / 6 - (7 * h[8] * b_formula[1]) / 6 + (31 * h[7] * b_formula[2]) / 48 - (37 * h[6] * b_formula[3]) / 144 + (23 * h[5] * b_formula[4]) / 288 -
                      (29 * h[4] * b_formula[5]) / 1440 + (73 * h[3] * b_formula[6]) / 17280 - (13 * h[2] * b_formula[7]) / 17280;
        Elem->Fy[6][5] = (-63 * h[10] * c_formula[0]) / 4 + (329 * h[9] * c_formula[1]) / 40 - (119 * h[8] * c_formula[2]) / 40 + (271 * h[7] * c_formula[3]) / 320 -
                      (197 * h[6] * c_formula[4]) / 960 + (17 * h[5] * c_formula[5]) / 384 - (83 * h[4] * c_formula[6]) / 9600 + (11 * h[3] * c_formula[7]) / 7200;
        Elem->Fy[6][6] = (-21 * h[11] * b_formula[0]) / 8 + (21 * h[10] * b_formula[1]) / 8 - (329 * h[9] * b_formula[2]) / 240 + (119 * h[8] * b_formula[3]) / 240 -
                      (271 * h[7] * b_formula[4]) / 1920 + (197 * h[6] * b_formula[5]) / 5760 - (17 * h[5] * b_formula[6]) / 2304 + (83 * h[4] * b_formula[7]) / 57600;
        Elem->Fy[6][7] = (165 * h[12] * c_formula[0]) / 4 - (339 * h[11] * c_formula[1]) / 16 + (119 * h[10] * c_formula[2]) / 16 - (193 * h[9] * c_formula[3]) / 96 +
                      (43 * h[8] * c_formula[4]) / 96 - (199 * h[7] * c_formula[5]) / 2304 + (1213 * h[6] * c_formula[6]) / 80640 - (11 * h[5] * c_formula[7]) / 4480;
        Elem->Fy[6][8] = (165 * h[13] * b_formula[0]) / 32 - (165 * h[12] * b_formula[1]) / 32 + (339 * h[11] * b_formula[2]) / 128 - (119 * h[10] * b_formula[3]) / 128 +
                      (193 * h[9] * b_formula[4]) / 768 - (43 * h[8] * b_formula[5]) / 768 + (199 * h[7] * b_formula[6]) / 18432 - (1213 * h[6] * b_formula[7]) / 645120;
        Elem->Fy[6][9] = (56485 * h[14] * c_formula[0]) / 384 - (55825 * h[13] * c_formula[1]) / 768 + (54505 * h[12] * c_formula[2]) / 2304 - (52435 * h[11] * c_formula[3]) / 9216 +
                      (9887 * h[10] * c_formula[4]) / 9216 - (27113 * h[9] * c_formula[5]) / 165888 + (23489 * h[8] * c_formula[6]) / 1161216 - (71 * h[7] * c_formula[7]) / 36288;
        Elem->Fy[6][10] = (11297 * h[15] * b_formula[0]) / 768 - (11297 * h[14] * b_formula[1]) / 768 + (11165 * h[13] * b_formula[2]) / 1536 - (10901 * h[12] * b_formula[3]) / 4608 +
                       (10487 * h[11] * b_formula[4]) / 18432 - (9887 * h[10] * b_formula[5]) / 92160 + (27113 * h[9] * b_formula[6]) / 1658880 - (23489 * h[8] * b_formula[7]) / 11612160;
        Elem->Fy[7][0] = b_formula[6] / 5040;
        Elem->Fy[7][1] = h[7] * c_formula[0] - (h[6] * c_formula[1]) / 2 + (h[5] * c_formula[2]) / 6 - (h[4] * c_formula[3]) / 24 + (h[3] * c_formula[4]) / 120 - (h[2] * c_formula[5]) / 720 +
                      (h[1] * c_formula[6]) / 5040 + c_formula[7] / 5040;
        Elem->Fy[7][2] = (h[8] * b_formula[0]) / 2 - (h[7] * b_formula[1]) / 2 + (h[6] * b_formula[2]) / 4 - (h[5] * b_formula[3]) / 12 + (h[4] * b_formula[4]) / 48 - (h[3] * b_formula[5]) / 240 +
                      (h[2] * b_formula[6]) / 1440 - (h[1] * b_formula[7]) / 10080;
        Elem->Fy[7][3] = -6 * h[9] * c_formula[0] + (13 * h[8] * c_formula[1]) / 4 - (5 * h[7] * c_formula[2]) / 4 + (3 * h[6] * c_formula[3]) / 8 - (11 * h[5] * c_formula[4]) / 120 +
                      (3 * h[4] * c_formula[5]) / 160 - (11 * h[3] * c_formula[6]) / 3360 + (h[2] * c_formula[7]) / 2016;
        Elem->Fy[7][4] = (-3 * h[10] * b_formula[0]) / 2 + (3 * h[9] * b_formula[1]) / 2 - (13 * h[8] * b_formula[2]) / 16 + (5 * h[7] * b_formula[3]) / 16 - (3 * h[6] * b_formula[4]) / 32 +
                      (11 * h[5] * b_formula[5]) / 480 - (3 * h[4] * b_formula[6]) / 640 + (11 * h[3] * b_formula[7]) / 13440;
        Elem->Fy[7][5] = (99 * h[11] * c_formula[0]) / 4 - (513 * h[10] * c_formula[1]) / 40 + (183 * h[9] * c_formula[2]) / 40 - (407 * h[8] * c_formula[3]) / 320 + (19 * h[7] * c_formula[4]) / 64 -
                      (39 * h[6] * c_formula[5]) / 640 + (757 * h[5] * c_formula[6]) / 67200 - (127 * h[4] * c_formula[7]) / 67200;
        Elem->Fy[7][6] = (33 * h[12] * b_formula[0]) / 8 - (33 * h[11] * b_formula[1]) / 8 + (171 * h[10] * b_formula[2]) / 80 - (61 * h[9] * b_formula[3]) / 80 + (407 * h[8] * b_formula[4]) / 1920 -
                      (19 * h[7] * b_formula[5]) / 384 + (13 * h[6] * b_formula[6]) / 1280 - (757 * h[5] * b_formula[7]) / 403200;
        Elem->Fy[7][7] = (4719 * h[13] * c_formula[0]) / 28 - (9339 * h[12] * c_formula[1]) / 112 + (3047 * h[11] * c_formula[2]) / 112 - (1471 * h[10] * c_formula[3]) / 224 +
                      (199 * h[9] * c_formula[4]) / 160 - (15331 * h[8] * c_formula[5]) / 80640 + (13213 * h[7] * c_formula[6]) / 564480 - (1229 * h[6] * c_formula[7]) / 564480;
        Elem->Fy[7][8] = (4719 * h[14] * b_formula[0]) / 224 - (4719 * h[13] * b_formula[1]) / 224 + (9339 * h[12] * b_formula[2]) / 896 - (3047 * h[11] * b_formula[3]) / 896 +
                      (1471 * h[10] * b_formula[4]) / 1792 - (199 * h[9] * b_formula[5]) / 1280 + (15331 * h[8] * b_formula[6]) / 645120 - (13213 * h[7] * b_formula[7]) / 4515840;
        Elem->Fy[7][9] = (6721 * h[15] * c_formula[0]) / 128 - (40755 * h[14] * c_formula[1]) / 1792 + (28171 * h[13] * c_formula[2]) / 5376 - (1375 * h[12] * c_formula[3]) / 3072 -
                      (4741 * h[11] * c_formula[4]) / 35840 + (14081 * h[10] * c_formula[5]) / 215040 - (691 * h[9] * c_formula[6]) / 43008 + (1819 * h[8] * c_formula[7]) / 645120;
        Elem->Fy[7][10] = (6721 * h[16] * b_formula[0]) / 1280 - (6721 * h[15] * b_formula[1]) / 1280 + (8151 * h[14] * b_formula[2]) / 3584 - (28171 * h[13] * b_formula[3]) / 53760 +
                       (275 * h[12] * b_formula[4]) / 6144 + (4741 * h[11] * b_formula[5]) / 358400 - (14081 * h[10] * b_formula[6]) / 2150400 + (691 * h[9] * b_formula[7]) / 430080;
        Elem->Fy[8][0] = b_formula[7] / 40320;
        Elem->Fy[8][1] = -(h[8] * c_formula[0]) + (h[7] * c_formula[1]) / 2 - (h[6] * c_formula[2]) / 6 + (h[5] * c_formula[3]) / 24 - (h[4] * c_formula[4]) / 120 + (h[3] * c_formula[5]) / 720 -
                      (h[2] * c_formula[6]) / 5040 + (h[1] * c_formula[7]) / 40320;
        Elem->Fy[8][2] = -(h[9] * b_formula[0]) / 2 + (h[8] * b_formula[1]) / 2 - (h[7] * b_formula[2]) / 4 + (h[6] * b_formula[3]) / 12 - (h[5] * b_formula[4]) / 48 + (h[4] * b_formula[5]) / 240 -
                      (h[3] * b_formula[6]) / 1440 + (h[2] * b_formula[7]) / 10080;
        Elem->Fy[8][3] = (15 * h[10] * c_formula[0]) / 2 - 4 * h[9] * c_formula[1] + (3 * h[8] * c_formula[2]) / 2 - (7 * h[7] * c_formula[3]) / 16 + (5 * h[6] * c_formula[4]) / 48 -
                      (h[5] * c_formula[5]) / 48 + (h[4] * c_formula[6]) / 280 - (43 * h[3] * c_formula[7]) / 80640;
        Elem->Fy[8][4] = (15 * h[11] * b_formula[0]) / 8 - (15 * h[10] * b_formula[1]) / 8 + h[9] * b_formula[2] - (3 * h[8] * b_formula[3]) / 8 + (7 * h[7] * b_formula[4]) / 64 -
                      (5 * h[6] * b_formula[5]) / 192 + (h[5] * b_formula[6]) / 192 - (h[4] * b_formula[7]) / 1120;
        Elem->Fy[8][5] = (-297 * h[12] * c_formula[0]) / 8 + (153 * h[11] * c_formula[1]) / 8 - (27 * h[10] * c_formula[2]) / 4 + (59 * h[9] * c_formula[3]) / 32 - (67 * h[8] * c_formula[4]) / 160 +
                      (53 * h[7] * c_formula[5]) / 640 - (197 * h[6] * c_formula[6]) / 13440 + (253 * h[5] * c_formula[7]) / 107520;
        Elem->Fy[8][6] = (-99 * h[13] * b_formula[0]) / 16 + (99 * h[12] * b_formula[1]) / 16 - (51 * h[11] * b_formula[2]) / 16 + (9 * h[10] * b_formula[3]) / 8 - (59 * h[9] * b_formula[4]) / 192 +
                      (67 * h[8] * b_formula[5]) / 960 - (53 * h[7] * b_formula[6]) / 3840 + (197 * h[6] * b_formula[7]) / 80640;
        Elem->Fy[8][7] = (-22737 * h[14] * c_formula[0]) / 112 + (2805 * h[13] * c_formula[1]) / 28 - (3641 * h[12] * c_formula[2]) / 112 + (3485 * h[11] * c_formula[3]) / 448 -
                      (3257 * h[10] * c_formula[4]) / 2240 + (4393 * h[9] * c_formula[5]) / 20160 - (367 * h[8] * c_formula[6]) / 14112 + (10291 * h[7] * c_formula[7]) / 4515840;
        Elem->Fy[8][8] = (-22737 * h[15] * b_formula[0]) / 896 + (22737 * h[14] * b_formula[1]) / 896 - (2805 * h[13] * b_formula[2]) / 224 + (3641 * h[12] * b_formula[3]) / 896 -
                      (3485 * h[11] * b_formula[4]) / 3584 + (3257 * h[10] * b_formula[5]) / 17920 - (4393 * h[9] * b_formula[6]) / 161280 + (367 * h[8] * b_formula[7]) / 112896;
        Elem->Fy[8][9] = (-24167 * h[16] * c_formula[0]) / 448 + (40755 * h[15] * c_formula[1]) / 1792 - (25597 * h[14] * c_formula[2]) / 5376 + (3355 * h[13] * c_formula[3]) / 21504 +
                      (8327 * h[12] * c_formula[4]) / 35840 - (57751 * h[11] * c_formula[5]) / 645120 + (18455 * h[10] * c_formula[6]) / 903168 - (13787 * h[9] * c_formula[7]) / 4014080;
        Elem->Fy[8][10] = (-24167 * h[17] * b_formula[0]) / 4480 + (24167 * h[16] * b_formula[1]) / 4480 - (8151 * h[15] * b_formula[2]) / 3584 +
                       (25597 * h[14] * b_formula[3]) / 53760 - (671 * h[13] * b_formula[4]) / 43008 - (8327 * h[12] * b_formula[5]) / 358400 + (57751 * h[11] * b_formula[6]) / 6451200 -
                       (3691 * h[10] * b_formula[7]) / 1806336;
        Elem->Fy[9][0] = 0;
        Elem->Fy[9][1] = h[9] * c_formula[0] - (h[8] * c_formula[1]) / 2 + (h[7] * c_formula[2]) / 6 - (h[6] * c_formula[3]) / 24 + (h[5] * c_formula[4]) / 120 - (h[4] * c_formula[5]) / 720 +
                      (h[3] * c_formula[6]) / 5040 - (h[2] * c_formula[7]) / 40320;
        Elem->Fy[9][2] = (h[10] * b_formula[0]) / 2 - (h[9] * b_formula[1]) / 2 + (h[8] * b_formula[2]) / 4 - (h[7] * b_formula[3]) / 12 + (h[6] * b_formula[4]) / 48 - (h[5] * b_formula[5]) / 240 +
                      (h[4] * b_formula[6]) / 1440 - (h[3] * b_formula[7]) / 10080;
        Elem->Fy[9][3] = (-55 * h[11] * c_formula[0]) / 6 + (29 * h[10] * c_formula[1]) / 6 - (16 * h[9] * c_formula[2]) / 9 + (73 * h[8] * c_formula[3]) / 144 - (17 * h[7] * c_formula[4]) / 144 +
                      (5 * h[6] * c_formula[5]) / 216 - (59 * h[5] * c_formula[6]) / 15120 + (139 * h[4] * c_formula[7]) / 241920;
        Elem->Fy[9][4] = (-55 * h[12] * b_formula[0]) / 24 + (55 * h[11] * b_formula[1]) / 24 - (29 * h[10] * b_formula[2]) / 24 + (4 * h[9] * b_formula[3]) / 9 - (73 * h[8] * b_formula[4]) / 576 +
                      (17 * h[7] * b_formula[5]) / 576 - (5 * h[6] * b_formula[6]) / 864 + (59 * h[5] * b_formula[7]) / 60480;
        Elem->Fy[9][5] = (-715 * h[13] * c_formula[0]) / 8 + 44 * h[12] * c_formula[1] - (341 * h[11] * c_formula[2]) / 24 + (323 * h[10] * c_formula[3]) / 96 - (59 * h[9] * c_formula[4]) / 96 +
                      (101 * h[8] * c_formula[5]) / 1152 - (379 * h[7] * c_formula[6]) / 40320 + (197 * h[6] * c_formula[7]) / 322560;
        Elem->Fy[9][6] = (-715 * h[14] * b_formula[0]) / 48 + (715 * h[13] * b_formula[1]) / 48 - (22 * h[12] * b_formula[2]) / 3 + (341 * h[11] * b_formula[3]) / 144 -
                      (323 * h[10] * b_formula[4]) / 576 + (59 * h[9] * b_formula[5]) / 576 - (101 * h[8] * b_formula[6]) / 6912 + (379 * h[7] * b_formula[7]) / 241920;
        Elem->Fy[9][7] = (-2145 * h[15] * c_formula[0]) / 112 + (715 * h[14] * c_formula[1]) / 112 - (1045 * h[12] * c_formula[3]) / 1344 + (473 * h[11] * c_formula[4]) / 1344 -
                      (193 * h[10] * c_formula[5]) / 2016 + (535 * h[9] * c_formula[6]) / 28224 - (2609 * h[8] * c_formula[7]) / 903168;
        Elem->Fy[9][8] = (-2145 * h[16] * b_formula[0]) / 896 + (2145 * h[15] * b_formula[1]) / 896 - (715 * h[14] * b_formula[2]) / 896 + (1045 * h[12] * b_formula[4]) / 10752 -
                      (473 * h[11] * b_formula[5]) / 10752 + (193 * h[10] * b_formula[6]) / 16128 - (535 * h[9] * b_formula[7]) / 225792;
        Elem->Fy[9][9] = (32175 * h[17] * c_formula[0]) / 448 - (9295 * h[16] * c_formula[1]) / 256 + (22165 * h[15] * c_formula[2]) / 1792 - (202345 * h[14] * c_formula[3]) / 64512 +
                      (39325 * h[13] * c_formula[4]) / 64512 - (103345 * h[12] * c_formula[5]) / 1161216 + (70477 * h[11] * c_formula[6]) / 8128512 - (12937 * h[10] * c_formula[7]) / 65028096;
        Elem->Fy[9][10] = (6435 * h[18] * b_formula[0]) / 896 - (6435 * h[17] * b_formula[1]) / 896 + (1859 * h[16] * b_formula[2]) / 512 - (4433 * h[15] * b_formula[3]) / 3584 +
                       (40469 * h[14] * b_formula[4]) / 129024 - (7865 * h[13] * b_formula[5]) / 129024 + (20669 * h[12] * b_formula[6]) / 2322432 - (70477 * h[11] * b_formula[7]) / 81285120;
        Elem->Fy[10][0] = 0;
        Elem->Fy[10][1] = -(h[10] * c_formula[0]) + (h[9] * c_formula[1]) / 2 - (h[8] * c_formula[2]) / 6 + (h[7] * c_formula[3]) / 24 - (h[6] * c_formula[4]) / 120 + (h[5] * c_formula[5]) / 720 -
                       (h[4] * c_formula[6]) / 5040 + (h[3] * c_formula[7]) / 40320;
        Elem->Fy[10][2] = -(h[11] * b_formula[0]) / 2 + (h[10] * b_formula[1]) / 2 - (h[9] * b_formula[2]) / 4 + (h[8] * b_formula[3]) / 12 - (h[7] * b_formula[4]) / 48 + (h[6] * b_formula[5]) / 240 -
                       (h[5] * b_formula[6]) / 1440 + (h[4] * b_formula[7]) / 10080;
        Elem->Fy[10][3] = 11 * h[12] * c_formula[0] - (23 * h[11] * c_formula[1]) / 4 + (25 * h[10] * c_formula[2]) / 12 - (7 * h[9] * c_formula[3]) / 12 + (2 * h[8] * c_formula[4]) / 15 -
                       (37 * h[7] * c_formula[5]) / 1440 + (43 * h[6] * c_formula[6]) / 10080 - (5 * h[5] * c_formula[7]) / 8064;
        Elem->Fy[10][4] = (11 * h[13] * b_formula[0]) / 4 - (11 * h[12] * b_formula[1]) / 4 + (23 * h[11] * b_formula[2]) / 16 - (25 * h[10] * b_formula[3]) / 48 + (7 * h[9] * b_formula[4]) / 48 -
                       (h[8] * b_formula[5]) / 30 + (37 * h[7] * b_formula[6]) / 5760 - (43 * h[6] * b_formula[7]) / 40320;
        Elem->Fy[10][5] = (3861 * h[14] * c_formula[0]) / 40 - (759 * h[13] * c_formula[1]) / 16 + (1221 * h[12] * c_formula[2]) / 80 - (115 * h[11] * c_formula[3]) / 32 +
                       (521 * h[10] * c_formula[4]) / 800 - (147 * h[9] * c_formula[5]) / 1600 + (13 * h[8] * c_formula[6]) / 1344 - (107 * h[7] * c_formula[7]) / 179200;
        Elem->Fy[10][6] = (1287 * h[15] * b_formula[0]) / 80 - (1287 * h[14] * b_formula[1]) / 80 + (253 * h[13] * b_formula[2]) / 32 - (407 * h[12] * b_formula[3]) / 160 +
                       (115 * h[11] * b_formula[4]) / 192 - (521 * h[10] * b_formula[5]) / 4800 + (49 * h[9] * b_formula[6]) / 3200 - (13 * h[8] * b_formula[7]) / 8064;
        Elem->Fy[10][7] = (1287 * h[16] * c_formula[0]) / 70 - (1287 * h[15] * c_formula[1]) / 224 - (429 * h[14] * c_formula[2]) / 1120 + (209 * h[13] * c_formula[3]) / 224 -
                       (1111 * h[12] * c_formula[4]) / 2800 + (21247 * h[11] * c_formula[5]) / 201600 - (5801 * h[10] * c_formula[6]) / 282240 + (17453 * h[9] * c_formula[7]) / 5644800;
        Elem->Fy[10][8] = (1287 * h[17] * b_formula[0]) / 560 - (1287 * h[16] * b_formula[1]) / 560 + (1287 * h[15] * b_formula[2]) / 1792 + (429 * h[14] * b_formula[3]) / 8960 -
                       (209 * h[13] * b_formula[4]) / 1792 + (1111 * h[12] * b_formula[5]) / 22400 - (21247 * h[11] * b_formula[6]) / 1612800 + (5801 * h[10] * b_formula[7]) / 2257920;
        Elem->Fy[10][9] = (-170599 * h[18] * c_formula[0]) / 2240 + (34463 * h[17] * c_formula[1]) / 896 - (175747 * h[16] * c_formula[2]) / 13440 +
                       (10153 * h[15] * c_formula[3]) / 3072 - (49049 * h[14] * c_formula[4]) / 76800 + (892309 * h[13] * c_formula[5]) / 9676800 - (118217 * h[12] * c_formula[6]) / 13547520 +
                       (72007 * h[11] * c_formula[7]) / 541900800;
        Elem->Fy[10][10] = (-170599 * h[19] * b_formula[0]) / 22400 + (170599 * h[18] * b_formula[1]) / 22400 - (34463 * h[17] * b_formula[2]) / 8960 +
                        (175747 * h[16] * b_formula[3]) / 134400 - (10153 * h[15] * b_formula[4]) / 30720 + (49049 * h[14] * b_formula[5]) / 768000 - (892309 * h[13] * b_formula[6]) / 96768000 +
                        (118217 * h[12] * b_formula[7]) / 135475200;

        Elem->Fx[1][1] = b_formula[1];
        Elem->Fx[0][2] = (h[1] * c_formula[0]) / 2 + c_formula[1] / 2;
        Elem->Fx[0][3] = (h[2] * b_formula[0]) / 6 - (h[1] * b_formula[1]) / 6 - b_formula[2] / 6;
        Elem->Fx[0][4] = -(h[3] * c_formula[0]) / 8 + (h[2] * c_formula[1]) / 8 - (h[1] * c_formula[2]) / 12 - c_formula[3] / 24;
        Elem->Fx[0][5] = -(h[4] * b_formula[0]) / 40 + (h[3] * b_formula[1]) / 40 - (h[2] * b_formula[2]) / 40 + (h[1] * b_formula[3]) / 60 + b_formula[4] / 120;
        Elem->Fx[0][6] = (h[5] * c_formula[0]) / 16 - (3 * h[4] * c_formula[1]) / 80 + (h[3] * c_formula[2]) / 60 - (h[2] * c_formula[3]) / 120 + (h[1] * c_formula[4]) / 240 + c_formula[5] / 720;
        Elem->Fx[0][7] = (h[6] * b_formula[0]) / 112 - (h[5] * b_formula[1]) / 112 + (3 * h[4] * b_formula[2]) / 560 - (h[3] * b_formula[3]) / 420 + (h[2] * b_formula[4]) / 840 -
                      (h[1] * b_formula[5]) / 1680 - b_formula[6] / 5040;
        Elem->Fx[0][8] = (-5 * h[7] * c_formula[0]) / 128 + (19 * h[6] * c_formula[1]) / 896 - (11 * h[5] * c_formula[2]) / 1344 + (h[4] * c_formula[3]) / 384 - (h[3] * c_formula[4]) / 1344 +
                      (h[2] * c_formula[5]) / 4032 - (h[1] * c_formula[6]) / 10080 - c_formula[7] / 40320;
        Elem->Fx[0][9] = (-5 * h[8] * b_formula[0]) / 1152 + (5 * h[7] * b_formula[1]) / 1152 - (19 * h[6] * b_formula[2]) / 8064 + (11 * h[5] * b_formula[3]) / 12096 -
                      (h[4] * b_formula[4]) / 3456 + (h[3] * b_formula[5]) / 12096 - (h[2] * b_formula[6]) / 36288 + (h[1] * b_formula[7]) / 90720;
        Elem->Fx[0][10] = (7 * h[9] * c_formula[0]) / 256 - (11 * h[8] * c_formula[1]) / 768 + (h[7] * c_formula[2]) / 192 - (h[6] * c_formula[3]) / 672 + (29 * h[5] * c_formula[4]) / 80640 -
                       (19 * h[4] * c_formula[5]) / 241920 + (h[3] * c_formula[6]) / 60480 - (h[2] * c_formula[7]) / 241920;
        Elem->Fx[1][2] = -(h[2] * c_formula[0]) + (h[1] * c_formula[1]) / 2 + c_formula[2] / 2;
        Elem->Fx[1][3] = -(h[3] * b_formula[0]) / 3 + (h[2] * b_formula[1]) / 3 - (h[1] * b_formula[2]) / 6 - b_formula[3] / 6;
        Elem->Fx[1][4] = (h[4] * c_formula[0]) / 2 - (3 * h[3] * c_formula[1]) / 8 + (5 * h[2] * c_formula[2]) / 24 - (h[1] * c_formula[3]) / 12 - c_formula[4] / 24;
        Elem->Fx[1][5] = (h[5] * b_formula[0]) / 10 - (h[4] * b_formula[1]) / 10 + (3 * h[3] * b_formula[2]) / 40 - (h[2] * b_formula[3]) / 24 + (h[1] * b_formula[4]) / 60 + b_formula[5] / 120;
        Elem->Fx[1][6] = (-3 * h[6] * c_formula[0]) / 8 + (17 * h[5] * c_formula[1]) / 80 - (7 * h[4] * c_formula[2]) / 80 + (h[3] * c_formula[3]) / 30 - (h[2] * c_formula[4]) / 80 +
                      (h[1] * c_formula[5]) / 240 + c_formula[6] / 720;
        Elem->Fx[1][7] = (-3 * h[7] * b_formula[0]) / 56 + (3 * h[6] * b_formula[1]) / 56 - (17 * h[5] * b_formula[2]) / 560 + (h[4] * b_formula[3]) / 80 - (h[3] * b_formula[4]) / 210 +
                      (h[2] * b_formula[5]) / 560 - (h[1] * b_formula[6]) / 1680 - b_formula[7] / 5040;
        Elem->Fx[1][8] = (5 * h[8] * c_formula[0]) / 16 - (149 * h[7] * c_formula[1]) / 896 + (167 * h[6] * c_formula[2]) / 2688 - (25 * h[5] * c_formula[3]) / 1344 +
                      (13 * h[4] * c_formula[4]) / 2688 - (5 * h[3] * c_formula[5]) / 4032 + (h[2] * c_formula[6]) / 2880 - (h[1] * c_formula[7]) / 10080;
        Elem->Fx[1][9] = (5 * h[9] * b_formula[0]) / 144 - (5 * h[8] * b_formula[1]) / 144 + (149 * h[7] * b_formula[2]) / 8064 - (167 * h[6] * b_formula[3]) / 24192 +
                      (25 * h[5] * b_formula[4]) / 12096 - (13 * h[4] * b_formula[5]) / 24192 + (5 * h[3] * b_formula[6]) / 36288 - (h[2] * b_formula[7]) / 25920;
        Elem->Fx[1][10] = (-35 * h[10] * c_formula[0]) / 128 + (109 * h[9] * c_formula[1]) / 768 - (13 * h[8] * c_formula[2]) / 256 + (19 * h[7] * c_formula[3]) / 1344 -
                       (53 * h[6] * c_formula[4]) / 16128 + (163 * h[5] * c_formula[5]) / 241920 - (31 * h[4] * c_formula[6]) / 241920 + (h[3] * c_formula[7]) / 40320;
        Elem->Fx[2][0] = -c_formula[1] / 2;
        Elem->Fx[2][1] = b_formula[2] / 2;
        Elem->Fx[2][2] = (3 * h[3] * c_formula[0]) / 2 - (3 * h[2] * c_formula[1]) / 4 + (h[1] * c_formula[2]) / 4 + c_formula[3] / 4;
        Elem->Fx[2][3] = (h[4] * b_formula[0]) / 2 - (h[3] * b_formula[1]) / 2 + (h[2] * b_formula[2]) / 4 - (h[1] * b_formula[3]) / 12 - b_formula[4] / 12;
        Elem->Fx[2][4] = (-5 * h[5] * c_formula[0]) / 4 + (13 * h[4] * c_formula[1]) / 16 - (19 * h[3] * c_formula[2]) / 48 + (7 * h[2] * c_formula[3]) / 48 - (h[1] * c_formula[4]) / 24 - c_formula[5] / 48;
        Elem->Fx[2][5] = -(h[6] * b_formula[0]) / 4 + (h[5] * b_formula[1]) / 4 - (13 * h[4] * b_formula[2]) / 80 + (19 * h[3] * b_formula[3]) / 240 - (7 * h[2] * b_formula[4]) / 240 +
                      (h[1] * b_formula[5]) / 120 + b_formula[6] / 240;
        Elem->Fx[2][6] = (21 * h[7] * c_formula[0]) / 16 - (23 * h[6] * c_formula[1]) / 32 + (9 * h[5] * c_formula[2]) / 32 - (3 * h[4] * c_formula[3]) / 32 + (7 * h[3] * c_formula[4]) / 240 -
                      (h[2] * c_formula[5]) / 120 + (h[1] * c_formula[6]) / 480 + c_formula[7] / 1440;
        Elem->Fx[2][7] = (3 * h[8] * b_formula[0]) / 16 - (3 * h[7] * b_formula[1]) / 16 + (23 * h[6] * b_formula[2]) / 224 - (9 * h[5] * b_formula[3]) / 224 + (3 * h[4] * b_formula[4]) / 224 -
                      (h[3] * b_formula[5]) / 240 + (h[2] * b_formula[6]) / 840 - (h[1] * b_formula[7]) / 3360;
        Elem->Fx[2][8] = (-45 * h[9] * c_formula[0]) / 32 + (189 * h[8] * c_formula[1]) / 256 - (69 * h[7] * c_formula[2]) / 256 + (139 * h[6] * c_formula[3]) / 1792 -
                      (17 * h[5] * c_formula[4]) / 896 + (23 * h[4] * c_formula[5]) / 5376 - (13 * h[3] * c_formula[6]) / 13440 + (h[2] * c_formula[7]) / 4480;
        Elem->Fx[2][9] = (-5 * h[10] * b_formula[0]) / 32 + (5 * h[9] * b_formula[1]) / 32 - (21 * h[8] * b_formula[2]) / 256 + (23 * h[7] * b_formula[3]) / 768 -
                      (139 * h[6] * b_formula[4]) / 16128 + (17 * h[5] * b_formula[5]) / 8064 - (23 * h[4] * b_formula[6]) / 48384 + (13 * h[3] * b_formula[7]) / 120960;
        Elem->Fx[2][10] = (385 * h[11] * c_formula[0]) / 256 - (397 * h[10] * c_formula[1]) / 512 + (421 * h[9] * c_formula[2]) / 1536 - (115 * h[8] * c_formula[3]) / 1536 +
                       (13 * h[7] * c_formula[4]) / 768 - (23 * h[6] * c_formula[5]) / 6912 + (41 * h[5] * c_formula[6]) / 69120 - (7 * h[4] * c_formula[7]) / 69120;
        Elem->Fx[3][0] = -c_formula[2] / 6;
        Elem->Fx[3][1] = b_formula[3] / 6;
        Elem->Fx[3][2] = -2 * h[4] * c_formula[0] + h[3] * c_formula[1] - (h[2] * c_formula[2]) / 3 + (h[1] * c_formula[3]) / 12 + c_formula[4] / 12;
        Elem->Fx[3][3] = (-2 * h[5] * b_formula[0]) / 3 + (2 * h[4] * b_formula[1]) / 3 - (h[3] * b_formula[2]) / 3 + (h[2] * b_formula[3]) / 9 - (h[1] * b_formula[4]) / 36 - b_formula[5] / 36;
        Elem->Fx[3][4] = (5 * h[6] * c_formula[0]) / 2 - (3 * h[5] * c_formula[1]) / 2 + (2 * h[4] * c_formula[2]) / 3 - (11 * h[3] * c_formula[3]) / 48 + (h[2] * c_formula[4]) / 16 - (h[1] * c_formula[5]) / 72 -
                      c_formula[6] / 144;
        Elem->Fx[3][5] = (h[7] * b_formula[0]) / 2 - (h[6] * b_formula[1]) / 2 + (3 * h[5] * b_formula[2]) / 10 - (2 * h[4] * b_formula[3]) / 15 + (11 * h[3] * b_formula[4]) / 240 -
                      (h[2] * b_formula[5]) / 80 + (h[1] * b_formula[6]) / 360 + b_formula[7] / 720;
        Elem->Fx[3][6] = (-7 * h[8] * c_formula[0]) / 2 + (15 * h[7] * c_formula[1]) / 8 - (17 * h[6] * c_formula[2]) / 24 + (7 * h[5] * c_formula[3]) / 32 - (29 * h[4] * c_formula[4]) / 480 +
                      (11 * h[3] * c_formula[5]) / 720 - (h[2] * c_formula[6]) / 288 + (h[1] * c_formula[7]) / 1440;
        Elem->Fx[3][7] = -(h[9] * b_formula[0]) / 2 + (h[8] * b_formula[1]) / 2 - (15 * h[7] * b_formula[2]) / 56 + (17 * h[6] * b_formula[3]) / 168 - (h[5] * b_formula[4]) / 32 +
                      (29 * h[4] * b_formula[5]) / 3360 - (11 * h[3] * b_formula[6]) / 5040 + (h[2] * b_formula[7]) / 2016;
        Elem->Fx[3][8] = (75 * h[10] * c_formula[0]) / 16 - (39 * h[9] * c_formula[1]) / 16 + (7 * h[8] * c_formula[2]) / 8 - (439 * h[7] * c_formula[3]) / 1792 + (103 * h[6] * c_formula[4]) / 1792 -
                      (97 * h[5] * c_formula[5]) / 8064 + (193 * h[4] * c_formula[6]) / 80640 - (19 * h[3] * c_formula[7]) / 40320;
        Elem->Fx[3][9] = (25 * h[11] * b_formula[0]) / 48 - (25 * h[10] * b_formula[1]) / 48 + (13 * h[9] * b_formula[2]) / 48 - (7 * h[8] * b_formula[3]) / 72 + (439 * h[7] * b_formula[4]) / 16128 -
                      (103 * h[6] * b_formula[5]) / 16128 + (97 * h[5] * b_formula[6]) / 72576 - (193 * h[4] * b_formula[7]) / 725760;
        Elem->Fx[3][10] = (-385 * h[12] * c_formula[0]) / 64 + (395 * h[11] * c_formula[1]) / 128 - (415 * h[10] * c_formula[2]) / 384 + (149 * h[9] * c_formula[3]) / 512 -
                       (33 * h[8] * c_formula[4]) / 512 + (85 * h[7] * c_formula[5]) / 6912 - (29 * h[6] * c_formula[6]) / 13824 + (23 * h[5] * c_formula[7]) / 69120;
        Elem->Fx[4][0] = -c_formula[3] / 24;
        Elem->Fx[4][1] = b_formula[4] / 24;
        Elem->Fx[4][2] = (5 * h[5] * c_formula[0]) / 2 - (5 * h[4] * c_formula[1]) / 4 + (5 * h[3] * c_formula[2]) / 12 - (5 * h[2] * c_formula[3]) / 48 + (h[1] * c_formula[4]) / 48 + c_formula[5] / 48;
        Elem->Fx[4][3] = (5 * h[6] * b_formula[0]) / 6 - (5 * h[5] * b_formula[1]) / 6 + (5 * h[4] * b_formula[2]) / 12 - (5 * h[3] * b_formula[3]) / 36 + (5 * h[2] * b_formula[4]) / 144 -
                      (h[1] * b_formula[5]) / 144 - b_formula[6] / 144;
        Elem->Fx[4][4] = (-35 * h[7] * c_formula[0]) / 8 + (5 * h[6] * c_formula[1]) / 2 - (25 * h[5] * c_formula[2]) / 24 + (65 * h[4] * c_formula[3]) / 192 - (17 * h[3] * c_formula[4]) / 192 +
                      (11 * h[2] * c_formula[5]) / 576 - (h[1] * c_formula[6]) / 288 - c_formula[7] / 576;
        Elem->Fx[4][5] = (-7 * h[8] * b_formula[0]) / 8 + (7 * h[7] * b_formula[1]) / 8 - (h[6] * b_formula[2]) / 2 + (5 * h[5] * b_formula[3]) / 24 - (13 * h[4] * b_formula[4]) / 192 +
                      (17 * h[3] * b_formula[5]) / 960 - (11 * h[2] * b_formula[6]) / 2880 + (h[1] * b_formula[7]) / 1440;
        Elem->Fx[4][6] = (63 * h[9] * c_formula[0]) / 8 - (133 * h[8] * c_formula[1]) / 32 + (49 * h[7] * c_formula[2]) / 32 - (173 * h[6] * c_formula[3]) / 384 + (221 * h[5] * c_formula[4]) / 1920 -
                      (17 * h[4] * c_formula[5]) / 640 + (h[3] * c_formula[6]) / 180 - (h[2] * c_formula[7]) / 960;
        Elem->Fx[4][7] = (9 * h[10] * b_formula[0]) / 8 - (9 * h[9] * b_formula[1]) / 8 + (19 * h[8] * b_formula[2]) / 32 - (7 * h[7] * b_formula[3]) / 32 + (173 * h[6] * b_formula[4]) / 2688 -
                      (221 * h[5] * b_formula[5]) / 13440 + (17 * h[4] * b_formula[6]) / 4480 - (h[3] * b_formula[7]) / 1260;
        Elem->Fx[4][8] = (-825 * h[11] * c_formula[0]) / 64 + (213 * h[10] * c_formula[1]) / 32 - (151 * h[9] * c_formula[2]) / 64 + (663 * h[8] * c_formula[3]) / 1024 -
                      (151 * h[7] * c_formula[4]) / 1024 + (271 * h[6] * c_formula[5]) / 9216 - (871 * h[5] * c_formula[6]) / 161280 + (307 * h[4] * c_formula[7]) / 322560;
        Elem->Fx[4][9] = (-275 * h[12] * b_formula[0]) / 192 + (275 * h[11] * b_formula[1]) / 192 - (71 * h[10] * b_formula[2]) / 96 + (151 * h[9] * b_formula[3]) / 576 -
                      (221 * h[8] * b_formula[4]) / 3072 + (151 * h[7] * b_formula[5]) / 9216 - (271 * h[6] * b_formula[6]) / 82944 + (871 * h[5] * b_formula[7]) / 1451520;
        Elem->Fx[4][10] = (-13299 * h[13] * c_formula[0]) / 256 + (13189 * h[12] * c_formula[1]) / 512 - (4323 * h[11] * c_formula[2]) / 512 + (4207 * h[10] * c_formula[3]) / 2048 -
                       (12109 * h[9] * c_formula[4]) / 30720 + (17051 * h[8] * c_formula[5]) / 276480 - (1927 * h[7] * c_formula[6]) / 241920 + (403 * h[6] * c_formula[7]) / 483840;
        Elem->Fx[5][0] = -c_formula[4] / 120;
        Elem->Fx[5][1] = b_formula[5] / 120;
        Elem->Fx[5][2] = -3 * h[6] * c_formula[0] + (3 * h[5] * c_formula[1]) / 2 - (h[4] * c_formula[2]) / 2 + (h[3] * c_formula[3]) / 8 - (h[2] * c_formula[4]) / 40 + (h[1] * c_formula[5]) / 240 +
                      c_formula[6] / 240;
        Elem->Fx[5][3] = -(h[7] * b_formula[0]) + h[6] * b_formula[1] - (h[5] * b_formula[2]) / 2 + (h[4] * b_formula[3]) / 6 - (h[3] * b_formula[4]) / 24 + (h[2] * b_formula[5]) / 120 -
                      (h[1] * b_formula[6]) / 720 - b_formula[7] / 720;
        Elem->Fx[5][4] = 7 * h[8] * c_formula[0] - (31 * h[7] * c_formula[1]) / 8 + (37 * h[6] * c_formula[2]) / 24 - (23 * h[5] * c_formula[3]) / 48 + (29 * h[4] * c_formula[4]) / 240 -
                      (73 * h[3] * c_formula[5]) / 2880 + (13 * h[2] * c_formula[6]) / 2880 - (h[1] * c_formula[7]) / 1440;
        Elem->Fx[5][5] = (7 * h[9] * b_formula[0]) / 5 - (7 * h[8] * b_formula[1]) / 5 + (31 * h[7] * b_formula[2]) / 40 - (37 * h[6] * b_formula[3]) / 120 + (23 * h[5] * b_formula[4]) / 240 -
                      (29 * h[4] * b_formula[5]) / 1200 + (73 * h[3] * b_formula[6]) / 14400 - (13 * h[2] * b_formula[7]) / 14400;
        Elem->Fx[5][6] = (-63 * h[10] * c_formula[0]) / 4 + (329 * h[9] * c_formula[1]) / 40 - (119 * h[8] * c_formula[2]) / 40 + (271 * h[7] * c_formula[3]) / 320 -
                      (197 * h[6] * c_formula[4]) / 960 + (17 * h[5] * c_formula[5]) / 384 - (83 * h[4] * c_formula[6]) / 9600 + (11 * h[3] * c_formula[7]) / 7200;
        Elem->Fx[5][7] = (-9 * h[11] * b_formula[0]) / 4 + (9 * h[10] * b_formula[1]) / 4 - (47 * h[9] * b_formula[2]) / 40 + (17 * h[8] * b_formula[3]) / 40 - (271 * h[7] * b_formula[4]) / 2240 +
                      (197 * h[6] * b_formula[5]) / 6720 - (17 * h[5] * b_formula[6]) / 2688 + (83 * h[4] * b_formula[7]) / 67200;
        Elem->Fx[5][8] = (495 * h[12] * c_formula[0]) / 16 - (1017 * h[11] * c_formula[1]) / 64 + (357 * h[10] * c_formula[2]) / 64 - (193 * h[9] * c_formula[3]) / 128 +
                      (43 * h[8] * c_formula[4]) / 128 - (199 * h[7] * c_formula[5]) / 3072 + (1213 * h[6] * c_formula[6]) / 107520 - (33 * h[5] * c_formula[7]) / 17920;
        Elem->Fx[5][9] = (55 * h[13] * b_formula[0]) / 16 - (55 * h[12] * b_formula[1]) / 16 + (113 * h[11] * b_formula[2]) / 64 - (119 * h[10] * b_formula[3]) / 192 +
                      (193 * h[9] * b_formula[4]) / 1152 - (43 * h[8] * b_formula[5]) / 1152 + (199 * h[7] * b_formula[6]) / 27648 - (1213 * h[6] * b_formula[7]) / 967680;
        Elem->Fx[5][10] = (11297 * h[14] * c_formula[0]) / 128 - (11165 * h[13] * c_formula[1]) / 256 + (10901 * h[12] * c_formula[2]) / 768 - (10487 * h[11] * c_formula[3]) / 3072 +
                       (9887 * h[10] * c_formula[4]) / 15360 - (27113 * h[9] * c_formula[5]) / 276480 + (23489 * h[8] * c_formula[6]) / 1935360 - (71 * h[7] * c_formula[7]) / 60480;
        Elem->Fx[6][0] = -c_formula[5] / 720;
        Elem->Fx[6][1] = b_formula[6] / 720;
        Elem->Fx[6][2] = (7 * h[7] * c_formula[0]) / 2 - (7 * h[6] * c_formula[1]) / 4 + (7 * h[5] * c_formula[2]) / 12 - (7 * h[4] * c_formula[3]) / 48 + (7 * h[3] * c_formula[4]) / 240 -
                      (7 * h[2] * c_formula[5]) / 1440 + (h[1] * c_formula[6]) / 1440 + c_formula[7] / 1440;
        Elem->Fx[6][3] = (7 * h[8] * b_formula[0]) / 6 - (7 * h[7] * b_formula[1]) / 6 + (7 * h[6] * b_formula[2]) / 12 - (7 * h[5] * b_formula[3]) / 36 + (7 * h[4] * b_formula[4]) / 144 -
                      (7 * h[3] * b_formula[5]) / 720 + (7 * h[2] * b_formula[6]) / 4320 - (h[1] * b_formula[7]) / 4320;
        Elem->Fx[6][4] = (-21 * h[9] * c_formula[0]) / 2 + (91 * h[8] * c_formula[1]) / 16 - (35 * h[7] * c_formula[2]) / 16 + (21 * h[6] * c_formula[3]) / 32 - (77 * h[5] * c_formula[4]) / 480 +
                      (21 * h[4] * c_formula[5]) / 640 - (11 * h[3] * c_formula[6]) / 1920 + (h[2] * c_formula[7]) / 1152;
        Elem->Fx[6][5] = (-21 * h[10] * b_formula[0]) / 10 + (21 * h[9] * b_formula[1]) / 10 - (91 * h[8] * b_formula[2]) / 80 + (7 * h[7] * b_formula[3]) / 16 - (21 * h[6] * b_formula[4]) / 160 +
                      (77 * h[5] * b_formula[5]) / 2400 - (21 * h[4] * b_formula[6]) / 3200 + (11 * h[3] * b_formula[7]) / 9600;
        Elem->Fx[6][6] = (231 * h[11] * c_formula[0]) / 8 - (1197 * h[10] * c_formula[1]) / 80 + (427 * h[9] * c_formula[2]) / 80 - (2849 * h[8] * c_formula[3]) / 1920 +
                      (133 * h[7] * c_formula[4]) / 384 - (91 * h[6] * c_formula[5]) / 1280 + (757 * h[5] * c_formula[6]) / 57600 - (127 * h[4] * c_formula[7]) / 57600;
        Elem->Fx[6][7] = (33 * h[12] * b_formula[0]) / 8 - (33 * h[11] * b_formula[1]) / 8 + (171 * h[10] * b_formula[2]) / 80 - (61 * h[9] * b_formula[3]) / 80 + (407 * h[8] * b_formula[4]) / 1920 -
                      (19 * h[7] * b_formula[5]) / 384 + (13 * h[6] * b_formula[6]) / 1280 - (757 * h[5] * b_formula[7]) / 403200;
        Elem->Fx[6][8] = (4719 * h[13] * c_formula[0]) / 32 - (9339 * h[12] * c_formula[1]) / 128 + (3047 * h[11] * c_formula[2]) / 128 - (1471 * h[10] * c_formula[3]) / 256 +
                      (1393 * h[9] * c_formula[4]) / 1280 - (15331 * h[8] * c_formula[5]) / 92160 + (13213 * h[7] * c_formula[6]) / 645120 - (1229 * h[6] * c_formula[7]) / 645120;
        Elem->Fx[6][9] = (1573 * h[14] * b_formula[0]) / 96 - (1573 * h[13] * b_formula[1]) / 96 + (3113 * h[12] * b_formula[2]) / 384 - (3047 * h[11] * b_formula[3]) / 1152 +
                      (1471 * h[10] * b_formula[4]) / 2304 - (1393 * h[9] * b_formula[5]) / 11520 + (15331 * h[8] * b_formula[6]) / 829440 - (13213 * h[7] * b_formula[7]) / 5806080;
        Elem->Fx[6][10] = (47047 * h[15] * c_formula[0]) / 1280 - (8151 * h[14] * c_formula[1]) / 512 + (28171 * h[13] * c_formula[2]) / 7680 - (1925 * h[12] * c_formula[3]) / 6144 -
                       (4741 * h[11] * c_formula[4]) / 51200 + (14081 * h[10] * c_formula[5]) / 307200 - (691 * h[9] * c_formula[6]) / 61440 + (1819 * h[8] * c_formula[7]) / 921600;
        Elem->Fx[7][0] = -c_formula[6] / 5040;
        Elem->Fx[7][1] = b_formula[7] / 5040;
        Elem->Fx[7][2] = -4 * h[8] * c_formula[0] + 2 * h[7] * c_formula[1] - (2 * h[6] * c_formula[2]) / 3 + (h[5] * c_formula[3]) / 6 - (h[4] * c_formula[4]) / 30 + (h[3] * c_formula[5]) / 180 -
                      (h[2] * c_formula[6]) / 1260 + (h[1] * c_formula[7]) / 10080;
        Elem->Fx[7][3] = (-4 * h[9] * b_formula[0]) / 3 + (4 * h[8] * b_formula[1]) / 3 - (2 * h[7] * b_formula[2]) / 3 + (2 * h[6] * b_formula[3]) / 9 - (h[5] * b_formula[4]) / 18 +
                      (h[4] * b_formula[5]) / 90 - (h[3] * b_formula[6]) / 540 + (h[2] * b_formula[7]) / 3780;
        Elem->Fx[7][4] = 15 * h[10] * c_formula[0] - 8 * h[9] * c_formula[1] + 3 * h[8] * c_formula[2] - (7 * h[7] * c_formula[3]) / 8 + (5 * h[6] * c_formula[4]) / 24 - (h[5] * c_formula[5]) / 24 +
                      (h[4] * c_formula[6]) / 140 - (43 * h[3] * c_formula[7]) / 40320;
        Elem->Fx[7][5] = 3 * h[11] * b_formula[0] - 3 * h[10] * b_formula[1] + (8 * h[9] * b_formula[2]) / 5 - (3 * h[8] * b_formula[3]) / 5 + (7 * h[7] * b_formula[4]) / 40 - (h[6] * b_formula[5]) / 24 +
                      (h[5] * b_formula[6]) / 120 - (h[4] * b_formula[7]) / 700;
        Elem->Fx[7][6] = (-99 * h[12] * c_formula[0]) / 2 + (51 * h[11] * c_formula[1]) / 2 - 9 * h[10] * c_formula[2] + (59 * h[9] * c_formula[3]) / 24 - (67 * h[8] * c_formula[4]) / 120 +
                      (53 * h[7] * c_formula[5]) / 480 - (197 * h[6] * c_formula[6]) / 10080 + (253 * h[5] * c_formula[7]) / 80640;
        Elem->Fx[7][7] = (-99 * h[13] * b_formula[0]) / 14 + (99 * h[12] * b_formula[1]) / 14 - (51 * h[11] * b_formula[2]) / 14 + (9 * h[10] * b_formula[3]) / 7 - (59 * h[9] * b_formula[4]) / 168 +
                      (67 * h[8] * b_formula[5]) / 840 - (53 * h[7] * b_formula[6]) / 3360 + (197 * h[6] * b_formula[7]) / 70560;
        Elem->Fx[7][8] = (-22737 * h[14] * c_formula[0]) / 112 + (2805 * h[13] * c_formula[1]) / 28 - (3641 * h[12] * c_formula[2]) / 112 + (3485 * h[11] * c_formula[3]) / 448 -
                      (3257 * h[10] * c_formula[4]) / 2240 + (4393 * h[9] * c_formula[5]) / 20160 - (367 * h[8] * c_formula[6]) / 14112 + (10291 * h[7] * c_formula[7]) / 4515840;
        Elem->Fx[7][9] = (-7579 * h[15] * b_formula[0]) / 336 + (7579 * h[14] * b_formula[1]) / 336 - (935 * h[13] * b_formula[2]) / 84 + (3641 * h[12] * b_formula[3]) / 1008 -
                      (3485 * h[11] * b_formula[4]) / 4032 + (3257 * h[10] * b_formula[5]) / 20160 - (4393 * h[9] * b_formula[6]) / 181440 + (367 * h[8] * b_formula[7]) / 127008;
        Elem->Fx[7][10] = (-24167 * h[16] * c_formula[0]) / 560 + (8151 * h[15] * c_formula[1]) / 448 - (25597 * h[14] * c_formula[2]) / 6720 + (671 * h[13] * c_formula[3]) / 5376 +
                       (8327 * h[12] * c_formula[4]) / 44800 - (57751 * h[11] * c_formula[5]) / 806400 + (3691 * h[10] * c_formula[6]) / 225792 - (13787 * h[9] * c_formula[7]) / 5017600;
        Elem->Fx[8][0] = -c_formula[7] / 40320;
        Elem->Fx[8][1] = 0;
        Elem->Fx[8][2] = (9 * h[9] * c_formula[0]) / 2 - (9 * h[8] * c_formula[1]) / 4 + (3 * h[7] * c_formula[2]) / 4 - (3 * h[6] * c_formula[3]) / 16 + (3 * h[5] * c_formula[4]) / 80 -
                      (h[4] * c_formula[5]) / 160 + (h[3] * c_formula[6]) / 1120 - (h[2] * c_formula[7]) / 8960;
        Elem->Fx[8][3] = (3 * h[10] * b_formula[0]) / 2 - (3 * h[9] * b_formula[1]) / 2 + (3 * h[8] * b_formula[2]) / 4 - (h[7] * b_formula[3]) / 4 + (h[6] * b_formula[4]) / 16 - (h[5] * b_formula[5]) / 80 +
                      (h[4] * b_formula[6]) / 480 - (h[3] * b_formula[7]) / 3360;
        Elem->Fx[8][4] = (-165 * h[11] * c_formula[0]) / 8 + (87 * h[10] * c_formula[1]) / 8 - 4 * h[9] * c_formula[2] + (73 * h[8] * c_formula[3]) / 64 - (17 * h[7] * c_formula[4]) / 64 +
                      (5 * h[6] * c_formula[5]) / 96 - (59 * h[5] * c_formula[6]) / 6720 + (139 * h[4] * c_formula[7]) / 107520;
        Elem->Fx[8][5] = (-33 * h[12] * b_formula[0]) / 8 + (33 * h[11] * b_formula[1]) / 8 - (87 * h[10] * b_formula[2]) / 40 + (4 * h[9] * b_formula[3]) / 5 - (73 * h[8] * b_formula[4]) / 320 +
                      (17 * h[7] * b_formula[5]) / 320 - (h[6] * b_formula[6]) / 96 + (59 * h[5] * b_formula[7]) / 33600;
        Elem->Fx[8][6] = (-2145 * h[13] * c_formula[0]) / 16 + 66 * h[12] * c_formula[1] - (341 * h[11] * c_formula[2]) / 16 + (323 * h[10] * c_formula[3]) / 64 - (59 * h[9] * c_formula[4]) / 64 +
                      (101 * h[8] * c_formula[5]) / 768 - (379 * h[7] * c_formula[6]) / 26880 + (197 * h[6] * c_formula[7]) / 215040;
        Elem->Fx[8][7] = (-2145 * h[14] * b_formula[0]) / 112 + (2145 * h[13] * b_formula[1]) / 112 - (66 * h[12] * b_formula[2]) / 7 + (341 * h[11] * b_formula[3]) / 112 -
                      (323 * h[10] * b_formula[4]) / 448 + (59 * h[9] * b_formula[5]) / 448 - (101 * h[8] * b_formula[6]) / 5376 + (379 * h[7] * b_formula[7]) / 188160;
        Elem->Fx[8][8] = (-19305 * h[15] * c_formula[0]) / 896 + (6435 * h[14] * c_formula[1]) / 896 - (3135 * h[12] * c_formula[3]) / 3584 + (1419 * h[11] * c_formula[4]) / 3584 -
                      (193 * h[10] * c_formula[5]) / 1792 + (535 * h[9] * c_formula[6]) / 25088 - (2609 * h[8] * c_formula[7]) / 802816;
        Elem->Fx[8][9] = (-2145 * h[16] * b_formula[0]) / 896 + (2145 * h[15] * b_formula[1]) / 896 - (715 * h[14] * b_formula[2]) / 896 + (1045 * h[12] * b_formula[4]) / 10752 -
                      (473 * h[11] * b_formula[5]) / 10752 + (193 * h[10] * b_formula[6]) / 16128 - (535 * h[9] * b_formula[7]) / 225792;
        Elem->Fx[8][10] = (57915 * h[17] * c_formula[0]) / 896 - (16731 * h[16] * c_formula[1]) / 512 + (39897 * h[15] * c_formula[2]) / 3584 - (40469 * h[14] * c_formula[3]) / 14336 +
                       (7865 * h[13] * c_formula[4]) / 14336 - (20669 * h[12] * c_formula[5]) / 258048 + (70477 * h[11] * c_formula[6]) / 9031680 - (12937 * h[10] * c_formula[7]) / 72253440;
        Elem->Fx[9][0] = 0;
        Elem->Fx[9][1] = 0;
        Elem->Fx[9][2] = -5 * h[10] * c_formula[0] + (5 * h[9] * c_formula[1]) / 2 - (5 * h[8] * c_formula[2]) / 6 + (5 * h[7] * c_formula[3]) / 24 - (h[6] * c_formula[4]) / 24 + (h[5] * c_formula[5]) / 144 -
                      (h[4] * c_formula[6]) / 1008 + (h[3] * c_formula[7]) / 8064;
        Elem->Fx[9][3] = (-5 * h[11] * b_formula[0]) / 3 + (5 * h[10] * b_formula[1]) / 3 - (5 * h[9] * b_formula[2]) / 6 + (5 * h[8] * b_formula[3]) / 18 - (5 * h[7] * b_formula[4]) / 72 +
                      (h[6] * b_formula[5]) / 72 - (h[5] * b_formula[6]) / 432 + (h[4] * b_formula[7]) / 3024;
        Elem->Fx[9][4] = (55 * h[12] * c_formula[0]) / 2 - (115 * h[11] * c_formula[1]) / 8 + (125 * h[10] * c_formula[2]) / 24 - (35 * h[9] * c_formula[3]) / 24 + (h[8] * c_formula[4]) / 3 -
                      (37 * h[7] * c_formula[5]) / 576 + (43 * h[6] * c_formula[6]) / 4032 - (25 * h[5] * c_formula[7]) / 16128;
        Elem->Fx[9][5] = (11 * h[13] * b_formula[0]) / 2 - (11 * h[12] * b_formula[1]) / 2 + (23 * h[11] * b_formula[2]) / 8 - (25 * h[10] * b_formula[3]) / 24 + (7 * h[9] * b_formula[4]) / 24 -
                      (h[8] * b_formula[5]) / 15 + (37 * h[7] * b_formula[6]) / 2880 - (43 * h[6] * b_formula[7]) / 20160;
        Elem->Fx[9][6] = (1287 * h[14] * c_formula[0]) / 8 - (1265 * h[13] * c_formula[1]) / 16 + (407 * h[12] * c_formula[2]) / 16 - (575 * h[11] * c_formula[3]) / 96 +
                      (521 * h[10] * c_formula[4]) / 480 - (49 * h[9] * c_formula[5]) / 320 + (65 * h[8] * c_formula[6]) / 4032 - (107 * h[7] * c_formula[7]) / 107520;
        Elem->Fx[9][7] = (1287 * h[15] * b_formula[0]) / 56 - (1287 * h[14] * b_formula[1]) / 56 + (1265 * h[13] * b_formula[2]) / 112 - (407 * h[12] * b_formula[3]) / 112 +
                      (575 * h[11] * b_formula[4]) / 672 - (521 * h[10] * b_formula[5]) / 3360 + (7 * h[9] * b_formula[6]) / 320 - (65 * h[8] * b_formula[7]) / 28224;
        Elem->Fx[9][8] = (1287 * h[16] * c_formula[0]) / 56 - (6435 * h[15] * c_formula[1]) / 896 - (429 * h[14] * c_formula[2]) / 896 + (1045 * h[13] * c_formula[3]) / 896 -
                      (1111 * h[12] * c_formula[4]) / 2240 + (21247 * h[11] * c_formula[5]) / 161280 - (5801 * h[10] * c_formula[6]) / 225792 + (17453 * h[9] * c_formula[7]) / 4515840;
        Elem->Fx[9][9] = (143 * h[17] * b_formula[0]) / 56 - (143 * h[16] * b_formula[1]) / 56 + (715 * h[15] * b_formula[2]) / 896 + (143 * h[14] * b_formula[3]) / 2688 -
                      (1045 * h[13] * b_formula[4]) / 8064 + (1111 * h[12] * b_formula[5]) / 20160 - (21247 * h[11] * b_formula[6]) / 1451520 + (5801 * h[10] * b_formula[7]) / 2032128;
        Elem->Fx[9][10] = (-170599 * h[18] * c_formula[0]) / 2240 + (34463 * h[17] * c_formula[1]) / 896 - (175747 * h[16] * c_formula[2]) / 13440 +
                       (10153 * h[15] * c_formula[3]) / 3072 - (49049 * h[14] * c_formula[4]) / 76800 + (892309 * h[13] * c_formula[5]) / 9676800 - (118217 * h[12] * c_formula[6]) / 13547520 +
                       (72007 * h[11] * c_formula[7]) / 541900800;
        Elem->Fx[10][0] = 0;
        Elem->Fx[10][1] = 0;
        Elem->Fx[10][2] = (11 * h[11] * c_formula[0]) / 2 - (11 * h[10] * c_formula[1]) / 4 + (11 * h[9] * c_formula[2]) / 12 - (11 * h[8] * c_formula[3]) / 48 + (11 * h[7] * c_formula[4]) / 240 -
                       (11 * h[6] * c_formula[5]) / 1440 + (11 * h[5] * c_formula[6]) / 10080 - (11 * h[4] * c_formula[7]) / 80640;
        Elem->Fx[10][3] = (11 * h[12] * b_formula[0]) / 6 - (11 * h[11] * b_formula[1]) / 6 + (11 * h[10] * b_formula[2]) / 12 - (11 * h[9] * b_formula[3]) / 36 + (11 * h[8] * b_formula[4]) / 144 -
                       (11 * h[7] * b_formula[5]) / 720 + (11 * h[6] * b_formula[6]) / 4320 - (11 * h[5] * b_formula[7]) / 30240;
        Elem->Fx[10][4] = (143 * h[13] * c_formula[0]) / 4 - (275 * h[12] * c_formula[1]) / 16 + (253 * h[11] * c_formula[2]) / 48 - (55 * h[10] * c_formula[3]) / 48 +
                       (11 * h[9] * c_formula[4]) / 60 - (121 * h[8] * c_formula[5]) / 5760 + (11 * h[7] * c_formula[6]) / 8064 + (11 * h[6] * c_formula[7]) / 161280;
        Elem->Fx[10][5] = (143 * h[14] * b_formula[0]) / 20 - (143 * h[13] * b_formula[1]) / 20 + (55 * h[12] * b_formula[2]) / 16 - (253 * h[11] * b_formula[3]) / 240 +
                       (11 * h[10] * b_formula[4]) / 48 - (11 * h[9] * b_formula[5]) / 300 + (121 * h[8] * b_formula[6]) / 28800 - (11 * h[7] * b_formula[7]) / 40320;
        Elem->Fx[10][6] = (-429 * h[15] * c_formula[0]) / 80 + (143 * h[14] * c_formula[1]) / 32 - (429 * h[13] * c_formula[2]) / 160 + (209 * h[12] * c_formula[3]) / 192 -
                       (1507 * h[11] * c_formula[4]) / 4800 + (649 * h[10] * c_formula[5]) / 9600 - (451 * h[9] * c_formula[6]) / 40320 + (1529 * h[8] * c_formula[7]) / 1075200;
        Elem->Fx[10][7] = (-429 * h[16] * b_formula[0]) / 560 + (429 * h[15] * b_formula[1]) / 560 - (143 * h[14] * b_formula[2]) / 224 + (429 * h[13] * b_formula[3]) / 1120 -
                       (209 * h[12] * b_formula[4]) / 1344 + (1507 * h[11] * b_formula[5]) / 33600 - (649 * h[10] * b_formula[6]) / 67200 + (451 * h[9] * b_formula[7]) / 282240;
        Elem->Fx[10][8] = (-21879 * h[17] * c_formula[0]) / 560 + (34749 * h[16] * c_formula[1]) / 1792 - (8151 * h[15] * c_formula[2]) / 1280 + (2717 * h[14] * c_formula[3]) / 1792 -
                       (5863 * h[13] * c_formula[4]) / 22400 + (46651 * h[12] * c_formula[5]) / 1612800 - (1133 * h[11] * c_formula[6]) / 2257920 - (26851 * h[10] * c_formula[7]) / 45158400;
        Elem->Fx[10][9] = (-2431 * h[18] * b_formula[0]) / 560 + (2431 * h[17] * b_formula[1]) / 560 - (3861 * h[16] * b_formula[2]) / 1792 + (2717 * h[15] * b_formula[3]) / 3840 -
                       (2717 * h[14] * b_formula[4]) / 16128 + (5863 * h[13] * b_formula[5]) / 201600 - (46651 * h[12] * b_formula[6]) / 14515200 + (1133 * h[11] * b_formula[7]) / 20321280;
        Elem->Fx[10][10] = (476333 * h[19] * c_formula[0]) / 22400 - (14443 * h[18] * c_formula[1]) / 1280 + (563849 * h[17] * c_formula[2]) / 134400 -
                        (260117 * h[16] * c_formula[3]) / 215040 + (1525381 * h[15] * c_formula[4]) / 5376000 - (5336903 * h[14] * c_formula[5]) / 96768000 +
                        (1191619 * h[13] * c_formula[6]) / 135475200 - (5947469 * h[12] * c_formula[7]) / 5419008000;
    }
}

static void compute_fields(const struct elem *Elem, double x, double y,
                           double *Fx, double *Fy)
{
    double xp[ELEGANT_MAX_EXPANSION_ORDER + 1];
    double yp[ELEGANT_MAX_EXPANSION_ORDER + 1];
    double sum_x = 0.0;
    double sum_y = 0.0;
    int i;
    int j;

    xp[0] = yp[0] = 1.0;
    for (i = 1; i <= Elem->ExpansionOrder; i++) {
        xp[i] = xp[i - 1] * x;
        yp[i] = yp[i - 1] * y;
    }

    for (i = 0; i <= Elem->ExpansionOrder; i++) {
        for (j = 0; j <= Elem->ExpansionOrder - i; j++) {
            sum_x += Elem->Fx[i][j] * xp[i] * yp[j];
            sum_y += Elem->Fy[i][j] * xp[i] * yp[j];
        }
    }
    *Fx = sum_x;
    *Fy = sum_y;
}

static int body_map(double *r6, const struct elem *Elem)
{
    const double *drift_frac;
    const double *kick_frac;
    double rho0 = Elem->Length / Elem->BendingAngle;
    double rho_actual = actual_bending_radius(Elem);
    double slice_length = Elem->Length / Elem->NumIntSteps;
    double distance = 0.0;
    int substeps;
    int i;
    int j;

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

    for (i = 0; i < Elem->NumIntSteps; i++) {
        for (j = 0; j < substeps; j++) {
            double dsh = slice_length * drift_frac[j];
            double f2 = (1.0 + r6[delta_]) * (1.0 + r6[delta_])
                      - r6[py_] * r6[py_];
            double f;
            double sin_phi;
            double phi;
            double sine;
            double cosine;
            double tangent;
            double cos_phi;
            double factor;

            if (f2 <= 0.0)
                return 0;
            f = sqrt(f2);
            sin_phi = r6[px_] / f;
            if (fabs(sin_phi) > 1.0)
                return 0;
            phi = asin(sin_phi);
            sine = sin(dsh / rho0 + phi);
            cosine = cos(dsh / rho0 + phi);
            if (cosine == 0.0)
                return 0;
            tangent = sine / cosine;
            cos_phi = cos(phi);

            r6[px_] = f * sine;
            factor = (rho0 + r6[x_]) * cos_phi / f
                   * (tangent - sin_phi / cos_phi);
            r6[y_] += r6[py_] * factor;
            distance += factor * (1.0 + r6[delta_]);
            f = cos_phi / cosine;
            r6[x_] = rho0 * (f - 1.0) + f * r6[x_];

            if (kick_frac[j] != 0.0) {
                double Fx;
                double Fy;
                double ds = slice_length * kick_frac[j];
                double kick_scale = ds * (1.0 + r6[x_] / rho0)
                                  / rho_actual;

                compute_fields(Elem, r6[x_], r6[y_], &Fx, &Fy);
                r6[px_] -= Fy * kick_scale;
                r6[py_] += Fx * kick_scale;
            }
        }
    }

    r6[ct_] += distance - Elem->Length;
    return 1;
}

static void elegant_csbend(double *r_in, const struct elem *Elem,
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

            apply_lindberg_fringe(r6, Elem, -1.0);

            if (!atIsNaN(r6[x_]) && !body_map(r6, Elem))
                r6[x_] = atGetNaN();

            if (!atIsNaN(r6[x_]) && Elem->FSECorrection == 1)
                r6[ct_] -= Elem->FSECorrectionPathError;

            if (!atIsNaN(r6[x_]))
                apply_lindberg_fringe(r6, Elem, 1.0);

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
    int requested_expansion;

    Elem->Length = atGetDouble(ElemData, "Length"); check_error();
    Elem->BendingAngle = atGetDouble(ElemData, "BendingAngle"); check_error();
    Elem->PolynomA = atGetDoubleArray(ElemData, "PolynomA"); check_error();
    Elem->PolynomB = atGetDoubleArray(ElemData, "PolynomB"); check_error();
    Elem->MaxOrder = atGetLong(ElemData, "MaxOrder"); check_error();
    Elem->NumIntSteps = atGetLong(ElemData, "NumIntSteps"); check_error();
    Elem->IntegrationOrder =
        atGetOptionalLong(ElemData, "IntegrationOrder", 4); check_error();
    Elem->Nonlinear =
        atGetOptionalLong(ElemData, "Nonlinear", 1); check_error();
    requested_expansion =
        atGetOptionalLong(ElemData, "ExpansionOrder", 0); check_error();
    Elem->Scaling =
        atGetOptionalDouble(ElemData, "FieldScaling", 1.0); check_error();
    Elem->EntranceAngle =
        atGetOptionalDouble(ElemData, "EntranceAngle", 0.0); check_error();
    Elem->ExitAngle =
        atGetOptionalDouble(ElemData, "ExitAngle", 0.0); check_error();
    Elem->FringeBendEntrance =
        atGetOptionalLong(ElemData, "FringeBendEntrance", 0); check_error();
    Elem->FringeBendExit =
        atGetOptionalLong(ElemData, "FringeBendExit", 0); check_error();
    Elem->FullGap =
        atGetOptionalDouble(ElemData, "FullGap", 0.0); check_error();
    Elem->FringeInt1 =
        atGetOptionalDouble(ElemData, "FringeInt1", 0.0); check_error();
    Elem->FringeInt2 =
        atGetOptionalDouble(ElemData, "FringeInt2", 0.0); check_error();
    Elem->H1 = atGetOptionalDouble(ElemData, "H1", 0.0); check_error();
    Elem->H2 = atGetOptionalDouble(ElemData, "H2", 0.0); check_error();
    Elem->FSE = atGetOptionalDouble(ElemData, "FSE", 0.0); check_error();
    Elem->FSEDipole =
        atGetOptionalDouble(ElemData, "FSE_DIPOLE", 0.0); check_error();
    Elem->FSEQuadrupole =
        atGetOptionalDouble(ElemData, "FSE_QUADRUPOLE", 0.0); check_error();
    Elem->FSECorrection =
        atGetOptionalLong(ElemData, "FSECorrection", 0); check_error();
    Elem->FSECorrection =
        atGetOptionalLong(ElemData, "FSE_CORRECTION",
                          Elem->FSECorrection); check_error();
    Elem->FSECorrectionValue =
        atGetOptionalDouble(ElemData, "FSECorrectionValue", 0.0);
    check_error();
    Elem->FSECorrectionValue =
        atGetOptionalDouble(ElemData, "FSE_CORRECTION_VALUE",
                            Elem->FSECorrectionValue); check_error();
    Elem->FSECorrectionPathError =
        atGetOptionalDouble(ElemData, "FSECorrectionPathError", 0.0);
    check_error();
    Elem->FSECorrectionPathError =
        atGetOptionalDouble(ElemData, "FSE_CORRECTION_PATH_ERROR",
                            Elem->FSECorrectionPathError); check_error();
    Elem->R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
    Elem->R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
    Elem->T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
    Elem->T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();
    Elem->EApertures =
        atGetOptionalDoubleArray(ElemData, "EApertures"); check_error();
    Elem->RApertures =
        atGetOptionalDoubleArray(ElemData, "RApertures"); check_error();

    if (Elem->Length == 0.0 || Elem->BendingAngle == 0.0)
        atError("ElegantCsbendPass: Length and BendingAngle must both be "
                "nonzero");
    if (Elem->NumIntSteps < 1)
        atError("ElegantCsbendPass: NumIntSteps must be positive");
    if (Elem->IntegrationOrder != 2 && Elem->IntegrationOrder != 4
            && Elem->IntegrationOrder != 6)
        atError("ElegantCsbendPass: IntegrationOrder must be 2, 4, or 6");
    if (Elem->MaxOrder < 0 || Elem->MaxOrder > ELEGANT_MAX_MULTIPOLE_ORDER)
        atError("ElegantCsbendPass: MaxOrder must be between 0 and 8");
    if (Elem->Scaling == 0.0)
        atError("ElegantCsbendPass: FieldScaling must be nonzero");

    Elem->ExpansionOrder = requested_expansion;
    if (Elem->ExpansionOrder < 0
            || Elem->ExpansionOrder > ELEGANT_MAX_EXPANSION_ORDER)
        atError("ElegantCsbendPass: ExpansionOrder must be between 0 and 10");
    check_error();

    compute_field_coefficients(Elem);
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

    elegant_csbend(r_in, Elem, num_particles);
    return Elem;
}

MODULE_DEF(ElegantCsbendPass)
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
        elegant_csbend(r_in, &Elem, num_particles);
    } else if (nrhs == 0) {
        plhs[0] = mxCreateCellMatrix(6, 1);
        mxSetCell(plhs[0], 0, mxCreateString("Length"));
        mxSetCell(plhs[0], 1, mxCreateString("BendingAngle"));
        mxSetCell(plhs[0], 2, mxCreateString("PolynomA"));
        mxSetCell(plhs[0], 3, mxCreateString("PolynomB"));
        mxSetCell(plhs[0], 4, mxCreateString("MaxOrder"));
        mxSetCell(plhs[0], 5, mxCreateString("NumIntSteps"));
        if (nlhs > 1) {
            plhs[1] = mxCreateCellMatrix(28, 1);
            mxSetCell(plhs[1], 0, mxCreateString("IntegrationOrder"));
            mxSetCell(plhs[1], 1, mxCreateString("Nonlinear"));
            mxSetCell(plhs[1], 2, mxCreateString("ExpansionOrder"));
            mxSetCell(plhs[1], 3, mxCreateString("FieldScaling"));
            mxSetCell(plhs[1], 4, mxCreateString("EntranceAngle"));
            mxSetCell(plhs[1], 5, mxCreateString("ExitAngle"));
            mxSetCell(plhs[1], 6, mxCreateString("FringeBendEntrance"));
            mxSetCell(plhs[1], 7, mxCreateString("FringeBendExit"));
            mxSetCell(plhs[1], 8, mxCreateString("FullGap"));
            mxSetCell(plhs[1], 9, mxCreateString("FringeInt1"));
            mxSetCell(plhs[1], 10, mxCreateString("FringeInt2"));
            mxSetCell(plhs[1], 11, mxCreateString("H1"));
            mxSetCell(plhs[1], 12, mxCreateString("H2"));
            mxSetCell(plhs[1], 13, mxCreateString("FSE"));
            mxSetCell(plhs[1], 14, mxCreateString("FSE_DIPOLE"));
            mxSetCell(plhs[1], 15, mxCreateString("FSE_QUADRUPOLE"));
            mxSetCell(plhs[1], 16, mxCreateString("FSECorrection"));
            mxSetCell(plhs[1], 17, mxCreateString("FSE_CORRECTION"));
            mxSetCell(plhs[1], 18, mxCreateString("FSECorrectionValue"));
            mxSetCell(plhs[1], 19, mxCreateString("FSE_CORRECTION_VALUE"));
            mxSetCell(plhs[1], 20, mxCreateString("FSECorrectionPathError"));
            mxSetCell(plhs[1], 21,
                      mxCreateString("FSE_CORRECTION_PATH_ERROR"));
            mxSetCell(plhs[1], 22, mxCreateString("R1"));
            mxSetCell(plhs[1], 23, mxCreateString("R2"));
            mxSetCell(plhs[1], 24, mxCreateString("T1"));
            mxSetCell(plhs[1], 25, mxCreateString("T2"));
            mxSetCell(plhs[1], 26, mxCreateString("EApertures"));
            mxSetCell(plhs[1], 27, mxCreateString("RApertures"));
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif
