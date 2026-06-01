#include <math.h>
#include "atlalib.c"
#define SQR(X) ((X)*(X))

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

static void yrot_propagate(double *r6, double c, double s, double pz, double p, double *bdiff)
{
  double dp1 = 1.0 + r6[delta_];
  double px = r6[px_];
  double py = r6[py_];
  double x_s_pz_p2 = r6[x_]*s/pz/p/p;
  double yrotmat[36];

  for (int m = 0; m < 36; m++)
    yrotmat[m] = 0.0;
  /* Set diagonal elements to 1	*/
  for (int m = 0; m < 6; m++)
    yrotmat[m * 7] = 1.0;

  yrotmat[0] = pz / p;           /* [0,0] */
  yrotmat[6] = x_s_pz_p2 * (SQR(px) + SQR(pz));  /* [0, 1] */
  yrotmat[18] = x_s_pz_p2 * px * py;          /* [0, 3] */
  yrotmat[24] = -x_s_pz_p2 * px * dp1;    /* [0, 4] */
  yrotmat[7] = c - s*px/pz;                       /* [1, 1] */
  yrotmat[19] = -s*py/pz;                        /* [1, 3] */
  yrotmat[25] = s*dp1/pz; /* [1, 4] */
  yrotmat[2] = s*py/p; /* [2, 0] */
  yrotmat[8] = x_s_pz_p2 * py * (c*px + s*pz); /* [2, 1] */
  yrotmat[20] = x_s_pz_p2 * (c*(SQR(pz)+SQR(py)) - s*px*pz);   /* [2, 3] */
  yrotmat[26] = -x_s_pz_p2 * c*py*dp1;  /* [2, 4] */
  yrotmat[5] = s*dp1/p; /* 5, 0] */
  yrotmat[11] = x_s_pz_p2 * dp1*(c*px + s*pz); /* [5, 1] */
  yrotmat[23] = x_s_pz_p2 * c*dp1*py; /* [5, 3] */
  yrotmat[29] = -x_s_pz_p2 * (c*(SQR(px)+SQR(py)) + s*px*pz); /* 5, 4] */

  ATsandwichmmt(yrotmat, bdiff);
}

static void Yrot(double *r6, double phi, double *bdiff)
{
    /* Forest 10.26, rotation in free space */

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

        if (bdiff) {
            yrot_propagate(r6, c, s, pz, p, bdiff);
        }
        r6[x_] = new_x;
        r6[px_] = new_px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }

static void bend_wedge_propagate(const double *r6, double irho, double fx, double fy, double fringecorr, double *bdiff)
{
    double p_norm = 1.0 / (1.0+r6[4]);

    for (int m = 0; m < 6; m++) {
    bdiff[1 + 6*m] += fx * bdiff[6*m];
    bdiff[3 + 6*m] -= fy * bdiff[2 + 6*m];
    }
    if (fringecorr != 0.0)
    for (int m = 0; m < 6; m++)
      bdiff[3 + 6*m] -= bdiff[4 + 6*m] * r6[2] * (irho*irho + fy*fy) * fringecorr * p_norm * p_norm / irho;

    for (int m = 0; m < 6; m++) {
    bdiff[m + 6*1] += fx * bdiff[m + 6*0];
    bdiff[m + 6*3] -= fy * bdiff[m + 6*2];
    }
    if (fringecorr != 0.0)
    for (int m = 0; m < 6; m++)
      bdiff[m + 6*3] -= bdiff[m + 6*4] * r6[2] * (irho*irho + fy*fy) * fringecorr * p_norm * p_norm / irho;
}

static void bend_linear_fringe(double* r6, double irho, double edge_angle,
        double gK, int method, double sign, double *bdiff)
{
    /*     method 0 no fringe field
     *     method 1 legacy version Brown First Order
     *     method 2 SOLEIL close to second order of Brown
     *     method 3 THOMX
     */
    double p_norm = 1.0 / (1.0+r6[4]);
    double fringecorr, fx, fy;

    /* Fringe field correction */
    if ((gK==0.0) || (method==0)) {
        fringecorr = 0.0;
    }
    else {
        register double sedge = sin(edge_angle);
        register double cedge = cos(edge_angle);
        fringecorr = irho * gK * (1.0 + sedge*sedge) / cedge;
    }

    /* Edge angle focusing */
    fx = irho * tan(edge_angle);
    if (method==1)
        fy = irho * tan(edge_angle - fringecorr*p_norm);
    else if (method==2)
        fy = irho * tan(edge_angle - fringecorr*p_norm) * p_norm;
    else if (method==3)
        fy = irho * tan(edge_angle - fringecorr + r6[1]*p_norm);
    else    /* fall back to legacy version */
        fy = irho * tan(edge_angle - fringecorr*p_norm);

    /*  Propagate B */
    if (bdiff) bend_wedge_propagate(r6, irho, fx, fy, fringecorr, bdiff);

  /*  Propagate particle */
    r6[px_] += r6[x_] * fx;
    r6[py_] -= r6[y_] * fy;
}

static void bend_wedge(double *r6, double rhoinv, double theta, double *bdiff)
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

        if (bdiff) {
            /* Propagation of the diffusion matrix temporarily computed as in bend_linear_fringe
               by using -theta */
            double fx = -rhoinv * tan(theta);
            double fy = fx;
            bend_wedge_propagate(r6, rhoinv, fx, fy, 0.0, bdiff);
        }

        r6[x_] = new_x;
        r6[px_] = new_px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }
    else {
        Yrot(r6, theta, bdiff);
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
