
static double Sec(double x)
{
    return 1.0 / cos(x);
}

#ifndef PXYZ
#define PXYZ
static double pxyz(double dp1, double px, double py)
{
  return sqrt(dp1*dp1 - px*px - py*py);
}
#endif /*PXYZ*/

static void Yrot(double *r6, double phi)
{
    /* Forest 10.26, rotation in free space */

    if (phi != 0.0) {
        double dp1 = 1.0 + r6[delta_];
        double c = cos(phi);
        double s = sin(phi);
        double pz = pxyz(dp1, r6[px_], r6[py_]);
        double p = c*pz - s*r6[px_];
        double px = s*pz + c*r6[px_];
        double x = r6[x_]*pz/p;
        double dy = r6[x_]*r6[py_]*s/p;
        double dct = dp1*r6[x_]*s/p;
        r6[x_] = x;
        r6[px_] = px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }
}

static void bend_fringe(double *r6, double irho, double gK)
{
    /* Forest 13.13, bend fringe in the hard-edge limit */

    double b0 = irho;

    double pz = pxyz(1.0+r6[delta_], r6[px_], r6[py_]);
    double px = r6[px_];
    double py = r6[py_];
    double d  = r6[delta_];
    double xp = px / pz;
    double yp = py / pz;

    double phi = -b0 * tan( b0 * gK * (1 + xp*xp*(2 + yp*yp))*pz - atan(xp / (1 + yp*yp)));

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
    double py2z2 = (py2 + pz2) * (py2 + pz2);
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
    double dxf = 0.5 * dpx * yf * yf;
    double dct = 0.5 * dd * yf * yf;
    double dpyf = phi * yf;

    r6[y_]  = yf;
    r6[x_]  += dxf;
    r6[py_] -= dpyf;
    r6[ct_] -= dct;
}

static void bend_edge(double *r6, double rhoinv, double theta)
{
    /* Forest 12.41, ideal wedge, map U(theta, rhoinv) */

    if (fabs(rhoinv) >= 1.e-6) {
        double dp1 = 1.0 + r6[4];
        double c = cos(theta);
        double s = sin(theta);
        double pz = pxyz(dp1, r6[px_], r6[py_]);
        double d2 = pxyz(dp1, 0.0, r6[py_]);
        double px = r6[px_]*c + (pz - rhoinv*r6[x_])*s;
        double dasin = asin(r6[px_]/d2) - asin(px/d2);
        double num = r6[x_]*(r6[px_]*sin(2.0*theta) + s*s*(2.0*pz - rhoinv*r6[x_]));
        double den = pxyz(dp1, px, r6[py_]) + pxyz(dp1, r6[px_], r6[py_])*c - r6[px_]*s;
        double x = r6[x_]*c + num/den;
        double dy = r6[py_]*theta/rhoinv + r6[py_]/rhoinv*dasin;
        double dct = dp1/rhoinv*(theta + dasin);

        r6[x_] = x;
        r6[px_] = px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }
}

static void bend_fringe_test(double *r6, double irho, double gK)
{
    /* Forest 13.13, bend fringe in the hard-edge limit */

    double b0 = irho;

    double pz = pxyz(1.0+r6[delta_], r6[px_], r6[py_]);
    double px = r6[px_];
    double py = r6[py_];
    double d  = r6[delta_];
    double xp = px / pz;
    double yp = py / pz;

    double psi = b0*xp/(1+yp*yp);

    /* these are the partial derivatives of psi with respect to px, py and delta
       total horror from Mathematica. This could benefit from some mini-TPSA */

    double px2 = px*px;
    double py2 = py*py;
    double pz2 = pz*pz;
    double pz4 = pz2*pz2;
    double denom = (pz2 + py2);

    double dpx = b0*(pz4 + pz2*denom - px2*py2)/pz/denom/denom;
    double dpy = -b0*yp*px/denom;
    double dd = b0*xp*(1.0+d)*(py2-pz2)/denom/denom;

    /* solve quadratic equation in yf (Forest fringe_part_I.pdf) */

    double yf = (2 * r6[y_]) / (1 + sqrt(1 - 2 * dpy * r6[y_]));
    double dxf = 0.5 * dpx * yf * yf;
    double dct = 0.5 * dd * yf * yf;
    double dpyf = psi * yf;
    // atPrintf("Fringe dx, dpy, dct: %g, %g, %g\n", dxf, dpyf, dct);


    r6[x_]  += dxf;
    r6[y_]  = yf;
    r6[py_] -= dpyf;
    r6[ct_] -= dct;
}
