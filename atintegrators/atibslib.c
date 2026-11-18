/*
 * Intra-beam scattering growth rates, shared by the IBSPass pass method and
 * the at.collective._ibs Python module.
 *
 * The optics are given as a (9, npoints) Fortran-ordered array, one column
 * per point: integration weight (normalised to 1), betax, betay, alphax,
 * alphay, Dx, D'x, Dy, D'y.
 *
 * Model: J.D. Bjorken, S.K. Mtingwa, Part. Accel. 13, 115 (1983), as modified
 * in MAD-X: F. Antoniou, F. Zimmermann, CERN-ATS-2012-066. Coulomb logarithm
 * computed as in MAD-X (twclog).
 */

#ifndef ATIBSLIB_C
#define ATIBSLIB_C

#include <math.h>
#include "atconstants.h"

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif

#define IBS_NOPTICS 9              /* weights, betx, bety, alfx, alfy, dx, dpx, dy, dpy */
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
    double bxbar = 0.0, bybar = 0.0, dxbar = 0.0, dybar = 0.0;
    int i;
    for (i = 0; i < npoints; i++) {
        const double *o = optics + IBS_NOPTICS*i;
        bxbar += o[0]*o[1];
        bybar += o[0]*o[2];
        dxbar += o[0]*o[5];
        dybar += o[0]*o[7];
    }
    double tempev = (gamma*ring->energy - ring->mass)*b->emitx/bxbar;
    double sigx = 100.0*sqrt(b->emitx*bxbar + dxbar*dxbar*b->sigma_e*b->sigma_e);
    double sigy = 100.0*sqrt(b->emity*bybar + dybar*dybar*b->sigma_e*b->sigma_e);
    double sigt = 100.0*b->bunch_length;
    double density = b->npart/(8.0*pow(M_PI, 1.5)*sigx*sigy*sigt);
    double debye = 743.4*sqrt(tempev/density)/fabs(ring->charge);
    double rmincl = 1.44e-7*ring->charge*ring->charge/tempev;
    double rminqm = __HBAR_C*1.0e5/(2.0*sqrt(1.0e-12*tempev*ring->mass));
    return log(fmin(sigx, debye)/fmax(rmincl, rminqm));
}

/* Ring averages of the integrals of Eq. (8) of the MAD-X note, including the
   factors in brackets */
static void ibs_bjorken_mtingwa(int npoints, const double *optics, double gamma,
                                const struct ibs_beam *b, double *integrals)
{
    double g2 = gamma*gamma;
    double s2 = 1.0/(b->sigma_e*b->sigma_e);
    double iex = 1.0/b->emitx, iey = 1.0/b->emity;
    int i, k;
    integrals[0] = integrals[1] = integrals[2] = 0.0;
    for (i = 0; i < npoints; i++) {
        const double *o = optics + IBS_NOPTICS*i;
        double betx = o[1], bety = o[2], alfx = o[3], alfy = o[4];
        double dx = o[5], dpx = o[6], dy = o[7], dpy = o[8];
        double gx = betx*iex, gy = bety*iey;
        double phix = dpx + alfx*dx/betx, phiy = dpy + alfy*dy/bety;
        double dx2 = dx*dx/betx, dy2 = dy*dy/bety;
        double hy0 = dy2 + bety*phiy*phiy;
        double hx = (dx2 + betx*phix*phix)*iex;
        double hy = hy0*iey;
        double ry = hy0/bety;
        double fx = gx*gx*phix*phix, fy = gy*gy*phiy*phiy;
        double hsum = hx + hy + s2;
        double dsum = g2*(dx2*iex + dy2*iey + s2);
        double a = g2*hsum + gx + gy;
        double bb = (gx+gy)*dsum + gx*gy*(g2*(phix*phix + phiy*phiy) + 1.0);
        double c = gx*gy*dsum;
        double ax = g2*hx*(2.0*g2*hsum - 2.0*gx - gy) - g2*gx*hy
                    + gx*(2.0*gx - gy - g2*s2) + 6.0*g2*fx;
        double bx = g2*hx*((gx+gy)*g2*hsum - g2*(fx+fy) + gx*(gx - 4.0*gy))
                    + gx*(g2*s2*(gx - 2.0*gy) + gx*gy*(1.0 + 6.0*g2*phix*phix)
                          + g2*(2.0*fy - fx))
                    + g2*gx*hy*(gx - 2.0*gy);
        double ay = gy*(-g2*(hx + 2.0*hy + gx*ry + s2)
                        + 2.0*g2*g2*ry*hsum - (gx - 2.0*gy) + 6.0*g2*gy*phiy*phiy);
        double by = gy*(g2*(gy - 2.0*gx)*(hx + s2) + g2*hy*(gy - 4.0*gx)
                        + gx*gy + g2*(2.0*fx - fy)
                        + g2*g2*ry*(gx+gy)*hsum - g2*g2*ry*(fx+fy)
                        + 6.0*g2*phiy*phiy*gx*gy);
        double az = g2*s2*(2.0*g2*hsum - gx - gy);
        double bz = g2*s2*((gx+gy)*g2*hsum - 2.0*gx*gy - g2*(fx+fy));
        double ix = 0.0, iy = 0.0, iz = 0.0;
        for (k = 0; k < IBS_NQUAD; k++) {
            double lam = ibs_lambda[k];
            double den = lam*(lam*(lam + a) + bb) + c;
            double f = ibs_weight[k]*sqrt(lam/den)/den;
            ix += f*(ax*lam + bx);
            iy += f*(ay*lam + by);
            iz += f*(az*lam + bz);
        }
        integrals[0] += o[0]*ix;
        integrals[1] += o[0]*iy;
        integrals[2] += o[0]*iz;
    }
}

/* Emittance growth rates [1/s]: horizontal, vertical and longitudinal
   (rate of sigma_delta^2) */
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
