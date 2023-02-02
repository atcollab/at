#include <math.h>

double Sec(double x)
{
    return 1.0 / cos(x);
}

void bend_fringe(double *r6, double irho, double gK)
{
    /* Forest 13.13, bend fringe in the hard-edge limit */

    double dpx, dpy, dd, b0, px, py, pz, d, phi, xp, yp, yf, xf, lf, pyf;

    b0 = irho;

    pz = get_pz(r6);
    px = r6[px_];
    py = r6[py_];
    d  = r6[delta_];
    xp = px / pz;
    yp = py / pz;

    phi = -b0 * tan( b0 * gK * (1 +    SQR(xp)*(2 + SQR(yp)))*pz - atan(xp / (1 + SQR(yp))));

    /* these are the partial derivatives of phi with respect to px, py and delta
       total horror from Mathematica. This could benefit from some mini-TPSA */

    dpx = -((b0*(pow(px,2)*pow(pz,4)*(pow(py,2) - pow(pz,2)) - pow(pz,6)*(pow(py,2) + pow(pz,2)) +
                 b0*gK*px*(pow(pz,2)*pow(pow(py,2) + pow(pz,2),2)*(2*pow(py,2) + 3*pow(pz,2)) + pow(px,4)*(3*pow(py,2)*pow(pz,2) + 2*pow(pz,4)) +
                           pow(px,2)*(3*pow(py,6) + 8*pow(py,4)*pow(pz,2) + 9*pow(py,2)*pow(pz,4) + 5*pow(pz,6))))*
             pow(Sec((b0*gK*(pow(pz,4) + pow(px,2)*(pow(py,2) + 2*pow(pz,2))))/pow(pz,3) - atan((px*pz)/(pow(py,2) + pow(pz,2)))),2))/
            (pow(pz,5)*(pow(py,4) + pow(px,2)*pow(pz,2) + 2*pow(py,2)*pow(pz,2) + pow(pz,4))));


    dpy = -((b0*py*(px*pow(pz,4)*(pow(py,2) + pow(pz,2)) + b0*gK*(-(pow(pz,4)*pow(pow(py,2) + pow(pz,2),2)) +
                                                                  pow(px,4)*(3*pow(py,2)*pow(pz,2) + 4*pow(pz,4)) +
                                                                  pow(px,2)*(3*pow(py,6) + 10*pow(py,4)*pow(pz,2) + 11*pow(py,2)*pow(pz,4) + 3*pow(pz,6))))*
             pow(Sec((b0*gK*(pow(pz,4) + pow(px,2)*(pow(py,2) + 2*pow(pz,2))))/pow(pz,3) - atan((px*pz)/(pow(py,2) + pow(pz,2)))),2))/
            (pow(pz,5)*(pow(py,4) + pow(px,2)*pow(pz,2) + 2*pow(py,2)*pow(pz,2) + pow(pz,4))));

    dd = (b0*(1 + d)*(px*pow(pz,4)*(pow(py,2) - pow(pz,2)) + b0*gK*
                      (-(pow(pz,4)*pow(pow(py,2) + pow(pz,2),2)) + pow(px,4)*(3*pow(py,2)*pow(pz,2) + 2*pow(pz,4)) +
                       pow(px,2)*(3*pow(py,6) + 8*pow(py,4)*pow(pz,2) + 7*pow(py,2)*pow(pz,4) + pow(pz,6))))*
          pow(Sec((b0*gK*(pow(pz,4) + pow(px,2)*(pow(py,2) + 2*pow(pz,2))))/pow(pz,3) - atan((px*pz)/(pow(py,2) + pow(pz,2)))),2))/
        (pow(pz,5)*(pow(py,4) + pow(px,2)*pow(pz,2) + 2*pow(py,2)*pow(pz,2) + pow(pz,4)));

    /* solve quadratic equation in yf (Forest fringe_part_I.pdf) */

    yf = (2 * r6[y_]) / (1 + sqrt(1 - 2 * dpy * r6[y_]));
    xf = r6[x_] + 0.5 * dpx * SQR(yf);
    lf = r6[ct_] - 0.5 * dd * SQR(yf);
    pyf = py - phi * yf;

    r6[y_]  = yf;
    r6[x_]  = xf;
    r6[py_] = pyf;
    r6[ct_] = lf;
}
