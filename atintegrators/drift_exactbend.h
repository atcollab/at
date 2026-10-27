#ifdef DIFFUSION
#error "drift_exactbend does not compute the diffusion matrix"
#endif

#include <math.h>
#include "atlalib.c"
/*
To speed up the integration loop, the path length is computed in absolute on
each step, and the reference total length is subtracted at the end of the loop.
For other uses, the relative path length is computed.
*/
#ifdef MAGNET_PASS
#define ABSOLUTE_PATH_LENGTH
#define FIX_LENGTH(length) r6[5] -= (length)
#endif

#ifndef PXYZ
#define PXYZ
static double pxyz(double dp1, double px, double py)
{
  return sqrt(dp1*dp1 - px*px - py*py);
}
#endif /*PXYZ*/

#define DRIFT(r6,length,irho,bdiff) drift(r6,length,irho)

static void drift(double *r6, double L, double irho)
{
    /* Forest 12.18, bend-kick split, map W(L,irho) */

    double dp1 = 1.0 + r6[delta_];
    double pz = pxyz(dp1, r6[px_], r6[py_]);
    if (fabs(irho) < 1.e-6) {
        double NormL = L / pz;
        r6[x_] += r6[px_] * NormL;
        r6[y_] += r6[py_] * NormL;
        r6[ct_] += NormL * dp1;   /* Absolute path length */
    }
    else {
        /* Rearranged to avoid cancellations between O(1) terms divided by irho.
           Half-angle: asin(px0/d) - asin(px/d) = 2 atan2(px0 - px, pz + pz1) */
        double rho = 1.0/irho;
        double px0 = r6[px_];
        double py = r6[py_];
        double sn = sin(irho*L);
        double cs = cos(irho*L);
        double shalf = sin(0.5*irho*L);
        double omc = 2.0*shalf*shalf;                   /* 1 - cos(irho*L) */
        double pzm1 = pz - 1.0;
        double q = pzm1 - r6[x_]*irho;
        double px = px0*cs + q*sn;
        double dpx = px0*omc - q*sn;                    /* px0 - px, accurate */
        double pz1 = pxyz(dp1, px, py);
        double spz = pz + pz1;
        double x = r6[x_]*(1.0-omc) + (px0*sn + pzm1*omc + dpx*(px0+px)/spz)*rho;
        double dasin = L + 2.0*atan2(dpx, spz)*rho;
        double dy = py*dasin;
        double dct = dp1*dasin;     /* Absolute path length */

        r6[x_] = x;
        r6[px_] = px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }
    #ifndef ABSOLUTE_PATH_LENGTH
    r6[ct_] -= L;
    #endif
}
