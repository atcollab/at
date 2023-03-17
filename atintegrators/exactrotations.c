#include <math.h>


static double denom(double dp1, double px, double py) {
  return sqrt(dp1*dp1 - px*px - py*py);
}

/* Forest 10.26, rotation in free space
   phi: y rotation
*/
static void Yrot(double phi, double *r6)
{
    double c, s;
    c = cos(phi);
    s = sin(phi);
    double r0[6] = {r6[0], r6[1], r6[2], r6[3], r6[4], r6[5]};
    double ps = get_pz(r6);
    double p = c*ps - s*r0[px_];
    r6[x_] = r0[x_]*ps/p;
    r6[px_] = s*ps + c*r0[px_];
    r6[y_] += r0[x_]*r0[py_]*s/p;
    r6[ct_] += (1.0+r0[delta_])*r0[x_]*s/p;
}

/* Forest 12.41, ideal wedge
   theta: angle
*/
static void Urot(double *r6, double rhoinv, double theta)
{
    if (theta != 0.0) {
        double r0[6] = {r6[0], r6[1], r6[2], r6[3], r6[4], r6[5]};
        double dp1 = 1.0 + r0[4];
        double c = cos(theta);
        double s = sin(theta);
        double ps = denom(dp1, r0[px_], r0[py_]);
        double dasin, num, den;
        double d2 = denom(dp1, 0.0, r0[py_]);

        r6[px_] = r0[px_]*c + (ps - rhoinv*r0[x_])*s;
        dasin = asin(r0[px_]/d2) - asin(r6[px_]/d2);
        num = r0[x_]*(r0[px_]*sin(2.0*theta) + s*s*(2.0*ps - rhoinv*r0[x_]));
        den = denom(dp1, r6[px_], r0[py_]) + denom(dp1, r0[px_], r0[py_])*c - r0[px_]*s;
        r6[x_] = r0[x_]*c + num/den;
        r6[y_] += r0[py_]*theta/rhoinv + r0[py_]/rhoinv*dasin;
        r6[ct_] += dp1/rhoinv*(theta + dasin);
    }
}
