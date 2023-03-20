#include <math.h>

double Sec(double x)
{
    return 1.0 / cos(x);
}

void bend_fringe(double *r6, double irho, double gK)
{
    /* Forest 13.13, bend fringe in the hard-edge limit */

    double b0 = irho;

    double pz = get_pz(r6);
    double px = r6[px_];
    double py = r6[py_];
    double d  = r6[delta_];
    double xp = px / pz;
    double yp = py / pz;

    double phi = -b0 * tan( b0 * gK * (1 + SQR(xp)*(2 + SQR(yp)))*pz - atan(xp / (1 + SQR(yp))));

    /* these are the partial derivatives of phi with respect to px, py and delta
       total horror from Mathematica. This could benefit from some mini-TPSA */

    double px2 = px*px;
    double px4 = px2*px2;
    double py2 = py*py;
    double py4 = py2*py2;
    double py6 = py4*py2;
    double pz2 = pz*pz;
    double pz3 = pz2*pz;
    double pz4 = pz2*pz2;
    double pz5 = pz4*pz;
    double pz6 = pz4*pz2;
    double py2z2 = SQR(py2 + pz2);
    double powsec = pow(Sec((b0*gK*(pz4 + px2*(py2 + 2*pz2)))/pz3 - atan((px*pz)/(py2 + pz2))),2);
    double denom = (pz5*(py4 + px2*pz2 + 2*py2*pz2 + pz4));

    double dpx = -(b0*(px2*pz4*(py2 - pz2) - pz6*(py2 + pz2) +
            b0*gK*px*(pz2*py2z2*(2*py2 + 3*pz2) + px4*(3*py2*pz2 + 2*pz4) +
                      px2*(3*py6 + 8*py4*pz2 + 9*py2*pz4 + 5*pz6)))*powsec)
           /denom;

    double dpy = -(b0*py*(px*pz4*(py2 + pz2) +
            b0*gK*(-(pz4*py2z2) + px4*(3*py2*pz2 + 4*pz4) +
                   px2*(3*py6 + 10*py4*pz2 + 11*py2*pz4 + 3*pz6)))*powsec)
           /denom;

    double dd = (b0*(1 + d)*(px*pz4*(py2 - pz2) + b0*gK*
                      (-(pz4*py2z2) + px4*(3*py2*pz2 + 2*pz4) +
                       px2*(3*py6 + 8*py4*pz2 + 7*py2*pz4 + pz6)))*powsec)
         /denom;

    /* solve quadratic equation in yf (Forest fringe_part_I.pdf) */

    double yf = (2 * r6[y_]) / (1 + sqrt(1 - 2 * dpy * r6[y_]));
    double xf = r6[x_] + 0.5 * dpx * SQR(yf);
    double lf = r6[ct_] - 0.5 * dd * SQR(yf);
    double pyf = py - phi * yf;

    r6[y_]  = yf;
    r6[x_]  = xf;
    r6[py_] = pyf;
    r6[ct_] = lf;
}

static double denom(double dp1, double px, double py) {
  return sqrt(dp1*dp1 - px*px - py*py);
}

/* Forest 12.41, ideal wedge: map U(theta,b1)
   theta: angle
*/
static void bend_edge(double *r6, double rhoinv, double theta)
{
/* Forest 12.41, ideal wedge, map U(theta, rhoinv) */
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
