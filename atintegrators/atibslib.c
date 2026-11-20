/*
 * Intra-beam scattering growth rates, shared by the IBSPass pass method and
 * the at.collective._ibs Python module.
 *
 * The optics are given as a (21, npoints) Fortran-ordered array, one column
 * per point: integration weight (normalised to 1), Dx, Dy, P1[6], P2[6],
 * P3[6]. Pk is the momentum block of J^T Rk J for the 6D mode k in the
 * coordinates (px, delta/gamma, py), the longitudinal mode being scaled so
 * that its emittance is sigma_delta^2.
 *
 * Model: J.D. Bjorken, S.K. Mtingwa, Part. Accel. 13, 115 (1983), as modified
 * in MAD-X: F. Antoniou, F. Zimmermann, CERN-ATS-2012-066. Coulomb logarithm
 * computed as in MAD-X (twclog). The local momentum distribution is built
 * from the 6D mode matrices, as in K. Kubo, K. Oide, PRSTAB 4, 124401 (2001).
 */

#ifndef ATIBSLIB_C
#define ATIBSLIB_C

#include <math.h>
#include "atconstants.h"

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif

#define IBS_NOPTICS 21             /* weight, dx, dy, P1[6], P2[6], P3[6] */
#define IBS_NDECADES 16            /* B&M integration over [1, 1e16] */
#define IBS_NGAUSS 8               /* Gauss-Legendre nodes per decade */
#define IBS_NQUAD (IBS_NDECADES*IBS_NGAUSS)

struct ibs_ring
{
    double energy;          /* Energy [eV] */
    double mass;            /* Rest energy [eV] */
    double charge;          /* Particle charge [e] */
    double beta;            /* Relativistic beta */
};

struct ibs_beam
{
    double npart;           /* Number of particles in the bunch */
    double emitx, emity, sigma_e, bunch_length;
};

static double ibs_lambda[IBS_NQUAD];
static double ibs_weight[IBS_NQUAD];
static int ibs_quadrature_ready = 0;

/* Nodes and weights of the Gauss-Legendre quadrature in ln(lambda),
   computed on first use */
static void ibs_quadrature(void)
{
    static const double x[IBS_NGAUSS] = {
        -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
    static const double w[IBS_NGAUSS] = {
        0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    int k;
    if (ibs_quadrature_ready) return;
    for (k = 0; k < IBS_NQUAD; k++) {
        ibs_lambda[k] = exp(M_LN10*(k/IBS_NGAUSS + 0.5*(x[k%IBS_NGAUSS] + 1.0)));
        ibs_weight[k] = 0.5*M_LN10*w[k%IBS_NGAUSS]*ibs_lambda[k];
    }
    ibs_quadrature_ready = 1;
}

/* Coulomb logarithm, as in MAD-X twclog */
static double ibs_coulomb_log(int npoints, const double *optics, double gamma,
                              const struct ibs_ring *ring, const struct ibs_beam *b)
{
    double b1x = 0.0, b1y = 0.0, b2x = 0.0, b2y = 0.0, dxbar = 0.0, dybar = 0.0;
    int i;
    for (i = 0; i < npoints; i++) {
        const double *o = optics + IBS_NOPTICS*i;
        b1x += o[0]*o[3];
        b1y += o[0]*o[8];
        b2x += o[0]*o[9];
        b2y += o[0]*o[14];
        dxbar += o[0]*o[1];
        dybar += o[0]*o[2];
    }
    double se2 = b->sigma_e*b->sigma_e;
    double tempev = (gamma*ring->energy - ring->mass)*b->emitx/b1x;
    double sigx = 100.0*sqrt(b->emitx*b1x + b->emity*b2x + dxbar*dxbar*se2);
    double sigy = 100.0*sqrt(b->emitx*b1y + b->emity*b2y + dybar*dybar*se2);
    double sigt = 100.0*b->bunch_length;
    double density = b->npart/(8.0*pow(M_PI, 1.5)*sigx*sigy*sigt);
    double debye = 743.4*sqrt(tempev/density)/fabs(ring->charge);
    double rmincl = 1.44e-7*ring->charge*ring->charge/tempev;
    double rminqm = __HBAR_C*1.0e5/(2.0*sqrt(1.0e-12*tempev*ring->mass));
    return log(fmin(sigx, debye)/fmax(rmincl, rminqm));
}

static double ibs_trace_sym(const double *a, const double *b)
{
    return a[0]*b[0] + a[3]*b[3] + a[5]*b[5] + 2.0*(a[1]*b[1] + a[2]*b[2] + a[4]*b[4]);
}

/* Ring averages of the integrals of Eq. (8) of the MAD-X note, including the
   factors in brackets */
static void ibs_bjorken_mtingwa(int npoints, const double *optics, double gamma,
                                const struct ibs_beam *b, double *integrals)
{
    double iemit[3] = {1.0/b->emitx, 1.0/b->emity, 1.0/(b->sigma_e*b->sigma_e)};
    int i, j, k;
    integrals[0] = integrals[1] = integrals[2] = 0.0;
    for (i = 0; i < npoints; i++) {
        const double *o = optics + IBS_NOPTICS*i;
        double lk[3][6], l[6] = {0.0}, adj[6], m[6];
        for (k = 0; k < 3; k++) {
            for (j = 0; j < 6; j++) {
                lk[k][j] = o[3+6*k+j]*iemit[k];
                l[j] += lk[k][j];
            }
        }
        adj[0] = l[3]*l[5] - l[4]*l[4];
        adj[1] = l[2]*l[4] - l[1]*l[5];
        adj[2] = l[1]*l[4] - l[2]*l[3];
        adj[3] = l[0]*l[5] - l[2]*l[2];
        adj[4] = l[1]*l[2] - l[0]*l[4];
        adj[5] = l[0]*l[3] - l[1]*l[1];
        double a = l[0] + l[3] + l[5];
        double bb = adj[0] + adj[3] + adj[5];
        double c = l[0]*adj[0] + l[1]*adj[1] + l[2]*adj[2];
        double i0 = 0.0, i1 = 0.0, diag;
        for (k = 0; k < IBS_NQUAD; k++) {
            double lam = ibs_lambda[k];
            double den = lam*(lam*(lam + a) + bb) + c;
            double f = ibs_weight[k]*sqrt(lam/den)/den;
            i0 += f;
            i1 += f*lam;
        }
        for (j = 0; j < 6; j++) m[j] = 3.0*(l[j]*i1 - adj[j]*i0);
        diag = bb*i0 - a*i1;
        m[0] += diag;
        m[3] += diag;
        m[5] += diag;
        for (k = 0; k < 3; k++) integrals[k] += o[0]*ibs_trace_sym(lk[k], m);
    }
}

/* Emittance growth rates [1/s] of the 3 modes */
static void ibs_growth_rates(int npoints, const double *optics,
                             const struct ibs_ring *ring,
                             const struct ibs_beam *b, double *rates)
{
    double gamma = ring->energy/ring->mass;
    double g2 = gamma*gamma;
    double beta3 = ring->beta*ring->beta*ring->beta;
    double r0 = __RE*ring->charge*ring->charge*__E0*1.0e9/ring->mass;
    double clog = ibs_coulomb_log(npoints, optics, gamma, ring, b);
    double cst = b->npart*r0*r0*C0*clog/(8.0*M_PI*beta3*g2*g2
                                         *b->emitx*b->emity*b->sigma_e*b->bunch_length);
    ibs_quadrature();
    ibs_bjorken_mtingwa(npoints, optics, gamma, b, rates);
    rates[0] *= cst;
    rates[1] *= cst;
    rates[2] *= cst;
}

#endif /*ATIBSLIB_C*/
