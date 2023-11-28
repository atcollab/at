#include <math.h>

#ifndef PXYZ
#define PXYZ
static double pxyz(double dp1, double px, double py)
{
  return sqrt(dp1*dp1 - px*px - py*py);
}
#endif /*PXYZ*/

static void exact_bend(double *r6, double irho, double L)
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
        double pzmx = pz-(1.0+r6[x_]*irho);
        double cs = cos(irho*L);
        double sn = sin(irho*L);
        double px = r6[px_]*cs + pzmx*sn;
        double d2 = pxyz(dp1, 0.0, r6[py_]);
        double dasin = L + (asin(r6[px_]/d2) - asin(px/d2))/irho;
        double x = (pxyz(dp1,px,r6[py_]) - pzmx*cs + r6[px_]*sn - 1.0)/irho;
        double dy = r6[py_]*dasin;
        double dct = dp1*dasin;     /* Absolute path length */

        r6[x_] = x;
        r6[px_] = px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }
}

static void exact_straight_bend(double *r6, double irho, double L)
{
    /* Forest 12.39, bend-kick split, map V(L,irho) */

    double dp1 = 1.0 + r6[delta_];
    double pz = pxyz(dp1, r6[px_], r6[py_]);
    if (fabs(irho) < 1.e-6) {
        double NormL = L / pz;
        r6[x_] += r6[px_] * NormL;
        r6[y_] += r6[py_] * NormL;
        r6[ct_] += NormL * dp1;   /* Absolute path length */
    }
    else {
        double px = r6[px_] - irho*L;
        double d2 = pxyz(dp1, 0.0, r6[py_]);
        double dasin = (asin(r6[px_]/d2) - asin(px/d2))/irho;
        double dx = (pxyz(dp1, px, r6[py_]) - pz)/irho;
        double dy = r6[py_]*dasin;
        double dct = dp1*dasin;     /* Absolute path length */

        r6[x_] += dx;
        r6[px_] = px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }
}
