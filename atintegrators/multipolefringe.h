#include "linearquadfringe.h"
#include "atlalib.c"

static void quad_fringe(double* r6, const double b2, double edge)
{
/*	Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam Dynamics..."by E. Forest */

  const double x = r6[x_];
  const double y = r6[y_];
  const double p_norm = 1.0 / (1.0 + r6[delta_]);
  const double u = edge * b2 / 12.0 * p_norm;
  const double x2 = x*x;
  const double y2 = y*y;
  const double xy = x*y;
  const double dx = u * (x2 + 3.0*y2) * x;
  const double dy = u * (y2 + 3.0*x2) * y;
/*
  const double dpx = 3.0 * u * (2.0*xy*r6[py_] - (x2+y2)*r6[px_]);
  const double dpy = 3.0 * u * (2.0*xy*r6[px_] - (x2+y2)*r6[py_]);

  r6[x_] += dx;
  r6[px_] += dpx;
  r6[y_] -= dy;
  r6[py_] -= dpy;
  r6[ct_] -= (dy*r6[3] - dx*r6[1]) * p_norm;
*/
  const double fxx = -edge * b2 / 4.0 * (x2+y2);
  const double fxy = -edge * b2 / 2.0 * xy;
  const double fyx = -fxy;
  const double fyy = -fxx;
  const double A = 1.0 - fxx * p_norm;
  const double B = -fyx * p_norm;
  const double C = -fxy * p_norm;
  const double D = 1.0 - fyy * p_norm;
  const double det = A*D - B*C;
  const double pxf = (D * r6[px_] - B * r6[py_]) / det;
  const double pyf = (A * r6[py_] - C * r6[px_]) / det;
  const double dct = (dx*pxf - dy*pyf) * p_norm;

  r6[x_] += dx;
  r6[px_] = pxf;
  r6[y_] -= dy;
  r6[py_] = pyf;
  r6[ct_] -= dct;
}

static void all_mult_fringe(double *r6,
                             const double *polya, const double *polyb, int max_order,
                             double edge)
{
  // Forest 13.29
  // note this is the sum over n of Forest 13.29
  // one for each multipole component

  const double x = r6[x_];
  const double y = r6[y_];
  const double p_norm = 1.0 / (1.0 + r6[delta_]);


  double FX = 0;
  double FY = 0;
  double FX_X = 0;
  double FX_Y = 0;
  double FY_X = 0;
  double FY_Y = 0;

  double RX = 1.0;
  double IX = 0.0;

  // invariant is (j is the index, i is the complex unit)
  // RX+IXi = (x + iy)^j
  for (int n = 0; n <= max_order; n++) {

    double DU = 0.0;
    double DV = 0.0;
    double U = 0.0;
    double V = 0.0;
    double B = polyb[n];
    double A = polya[n];

    double j = n + 1.0;
    double f1 = -edge / 4.0 / (j + 1.0);
    double nf = 1.0 * (j + 2.0) / j;

    double DRX = RX;
    double DIX = IX;

    // complex muls

    RX = DRX * x - DIX * y;
    IX = DRX * y + DIX * x;

    if (n > 0) {  // Skip field order 0
      U = f1 * (B * RX - A * IX);
      V = f1 * (B * IX + A * RX);
      DU = f1 * (B * DRX - A * DIX);
      DV = f1 * (B * DIX + A * DRX);
    }

    const double DUX = j * DU;
    const double DVX = j * DV;
    const double DUY = -j * DV;
    const double DVY = j * DU;


    FX += U * x + nf * V * y;
    FY += U * y - nf * V * x;

    FX_X += DUX * x + U + nf * y * DVX;
    FX_Y += DUY * x + nf * (V + y * DVY);
    FY_X += DUX * y - nf * (V + x * DVX);
    FY_Y += DUY * y + U - nf * x * DVY;
  }


  // solve 2x2 matrix equation

  double A = 1.0 - FX_X * p_norm;
  double B = -FY_X * p_norm;
  double C = -FX_Y * p_norm;
  double D = 1.0 - FY_Y * p_norm;
  double det = A*D - B*C;

  double pxf = (D * r6[px_] - B * r6[py_]) / det;
  double pyf = (A * r6[py_] - C * r6[px_]) / det;

  r6[x_] -= FX * p_norm;
  r6[y_] -= FY * p_norm;
  r6[px_] = pxf;
  r6[py_] = pyf;
  r6[ct_] = r6[ct_] - (pxf * FX + pyf * FY) * p_norm * p_norm;
}

static void multipole_fringe(double *r6, int method, double B1,
                             const double *polya, const double *polyb, int max_order,
                             double *fringeIntM0, double *fringeIntP0,
                             double edge)
{
    switch (method) {
        case 0: break;
        case 1: quad_fringe(r6, B1, edge);
                break;
        case 2: if (fringeIntM0 && fringeIntP0) {
                    if (edge > 0)
                        linearQuadFringeElegantEntrance(r6, B1, fringeIntM0, fringeIntP0);
                    else
                        linearQuadFringeElegantExit(r6, B1, fringeIntM0, fringeIntP0);
                }
                else
                    quad_fringe(r6, B1, edge);
                break;
        case 3: all_mult_fringe(r6, polya, polyb, max_order, edge);
                break;
    }
}