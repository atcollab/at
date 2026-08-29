
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
        double x = r6[x_];
        double px = r6[px_];
        double py = r6[py_];
        double c = cos(phi);
        double s = sin(phi);
        double pz = pxyz(dp1, px, py);
        double p = c*pz - s*px;
        double new_px = s*pz + c*px;
        double new_x = x*pz/p;
        double dy = x*py*s/p;
        double dct = dp1*x*s/p;
        r6[x_] = new_x;
        r6[px_] = new_px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }
}

static void bend_fringe_orig(double *r6, double irho, double gK)
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
        double x = r6[x_];
        double px = r6[px_];
        double py = r6[py_];
        double pz = pxyz(dp1, px, py);
        double d2 = pxyz(dp1, 0.0, py);
        double new_px = px*c + (pz - rhoinv*x)*s;
        double dasin = asin(px/d2) - asin(new_px/d2);
        double num = x*s*(2.0*px*c + s*(2.0*pz - rhoinv*x));
        double den = pxyz(dp1, new_px, py) + pz*c - px*s;
        double new_x = x*c + num/den;
        double dy = py*theta/rhoinv + py/rhoinv*dasin;
        double dct = dp1/rhoinv*(theta + dasin);

        r6[x_] = new_x;
        r6[px_] = new_px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }
    else {
        Yrot(r6, theta);
    }
}

static void quad_wedge(double *r6, double k1_theta)
{
    double x = r6[x_];
    double y = r6[y_];
    double dpx = k1_theta * (x*x - 0.5*y*y);
    double dpy = k1_theta * x * y;
    r6[px_] -= dpx;
    r6[py_] += dpy;
}

void bend_fringe(double *r6, double irho, double gK)
{
    double factor = 0.0;
    if (fabs(gK) > 1.0e-6) {
        factor = SQR(irho) / 9.0 / gK;
    }
    const double irho_g_fint = irho * gK;
    const double y = r6[y_];
    const double px = r6[px_];
    const double py = r6[py_];
    const double dp1 = r6[delta_] + 1.0;

    const double pz2 = SQR(dp1) - SQR(px) - SQR(py);
    const double pz = sqrt(pz2);
    const double xp = px / pz;
    const double yp = py / pz;

    const double dpz_dpx = -xp;
    const double dpz_dpy = -yp;
    const double dpz_ddelta = dp1 / pz;

    const double dxp_dpx =    -px/pz2 * dpz_dpx     + 1/pz;
    const double dxp_dpy =    -px/pz2 * dpz_dpy;
    const double dxp_ddelta = -px/pz2 * dpz_ddelta;

    const double dyp_dpx =    -py/pz2 * dpz_dpx;
    const double dyp_dpy =    -py/pz2 * dpz_dpy     + 1/pz;
    const double dyp_ddelta = -py/pz2 * dpz_ddelta;

    const double phi0 = xp / (1.0 + SQR(yp));
    const double dphi0_dxp = 1.0 / (1.0 + SQR(yp));
    const double dphi0_dyp = -2 * xp * yp / SQR(1.0 + SQR(yp));

    const double dphi0_dpx = dphi0_dxp * dxp_dpx + dphi0_dyp * dyp_dpx;
    const double dphi0_dpy = dphi0_dxp * dxp_dpy + dphi0_dyp * dyp_dpy;
    const double dphi0_ddelta = dphi0_dxp * dxp_ddelta + dphi0_dyp * dyp_ddelta;

    const double phi1 = 1.0 + 2.0 * SQR(xp) + SQR(xp) * SQR(yp);
    const double dphi1_dxp = 4.0 * xp + 2.0 * SQR(yp) * xp;
    const double dphi1_dyp = 2.0 * SQR(xp) * yp;

    const double dphi1_dpx = dphi1_dxp * dxp_dpx + dphi1_dyp * dyp_dpx;
    const double dphi1_dpy = dphi1_dxp * dxp_dpy + dphi1_dyp * dyp_dpy;
    const double dphi1_ddelta = dphi1_dxp * dxp_ddelta + dphi1_dyp * dyp_ddelta;

    const double phi2 = atan(phi0) - irho_g_fint * pz * phi1;
    const double dphi2_dpx = dphi0_dpx / (1.0 + SQR(phi0))
                 - irho_g_fint * (pz * dphi1_dpx + phi1 * dpz_dpx);
    const double dphi2_dpy = dphi0_dpy / (1.0 + SQR(phi0))
                 - irho_g_fint * (pz * dphi1_dpy + phi1 * dpz_dpy);
    const double dphi2_ddelta = dphi0_ddelta / (1.0 + SQR(phi0))
                 - irho_g_fint * (pz * dphi1_ddelta + phi1 * dpz_ddelta);

    const double Phi0 = irho * tan(phi2);
    const double cphi2 = cos(phi2);
    const double irho_c2 = irho / cphi2 / cphi2;
    const double dPhi0_dpx = irho_c2 * dphi2_dpx;
    const double dPhi0_dpy = irho_c2 * dphi2_dpy;
    const double dPhi0_ddelta = irho_c2 * dphi2_ddelta;

    const double new_y = 2.0 * y / (1.0 + sqrt(1.0 - 2.0 * dPhi0_dpy * y));
    const double delta_x = dPhi0_dpx * SQR(new_y) / 2;
    const double delta_py = -Phi0 * new_y - factor / dp1 * SQR(new_y) * new_y;
    const double delta_l = -dPhi0_ddelta * SQR(new_y) / 2.0;

    r6[x_]  += delta_x;
    r6[y_]  = new_y;
    r6[py_] += delta_py;
    r6[ct_] += delta_l;
}
