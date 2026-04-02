
static void multipole_fringe(double *r6, double L,
                             double *polya, double *polyb, int max_order,
                             double edge, int skip_b0)
{
  // PTC multipole_fringer
  // Forest 13.29
  // not re-derived and checked
  // note this is the sum over n of Forest 13.29
  // one for each multipole component

  double U, V, DU, DV, DUX, DVX, DUY, DVY, FX, FY, FX_X, FX_Y, FY_X, FY_Y,
      RX, IX, DRX, DIX;

  FX = 0;
  FY = 0;
  FX_X = 0;
  FX_Y = 0;
  FY_X = 0;
  FY_Y = 0;

  RX = 1.0;
  IX = 0.0;

  // invariant is (j is the index, i is the complex unit)
  // RX+IXi = (x + iy)^j
  for (int n = 0; n <= max_order; n++) {

    double B = polyb[n];
    double A = polya[n];

    int j = n + 1;

    DRX = RX;
    DIX = IX;

    // complex muls

    RX = DRX * r6[x_] - DIX * r6[y_];
    IX = DRX * r6[y_] + DIX * r6[x_];

    if (n == 0 && skip_b0) {
      U  =         - A * IX;
      V  =         + A * RX;
      DU =         - A * DIX;
      DV =         + A * DRX;
    }
    else {
      U = B * RX - A * IX;
      V = B * IX + A * RX;
      DU = B * DRX - A * DIX;
      DV = B * DIX + A * DRX;
    }
    double f1 = -edge / 4.0 / (j + 1);

    U = U * f1;
    V = V * f1;
    DU = DU * f1;
    DV = DV * f1;

    DUX = j * DU;
    DVX = j * DV;
    DUY = -j * DV;
    DVY = j * DU;

    double nf = 1.0 * (j + 2) / j;

    FX += U * r6[x_] + nf * V * r6[y_];
    FY += U * r6[y_] - nf * V * r6[x_];

    FX_X += DUX * r6[x_] + U + nf * r6[y_] * DVX;
    FX_Y += DUY * r6[x_] + nf * V + nf * r6[y_] * DVY;

    FY_X += DUX * r6[y_] - nf * V - nf * r6[x_] * DVX;
    FY_Y += DUY * r6[y_] + U - nf * r6[x_] * DVY;
  }

  double DEL = 1.0 / (1 + r6[delta_]);

  // solve 2x2 matrix equation

  double A = 1 - FX_X * DEL;
  double B = -FY_X * DEL;
  double D = 1 - FY_Y * DEL;
  double C = -FX_Y * DEL;

  r6[x_] = r6[x_] - FX * DEL;
  r6[y_] = r6[y_] - FY * DEL;

  double pxf = (D * r6[px_] - B * r6[py_]) / (A * D - B * C);
  double pyf = (A * r6[py_] - C * r6[px_]) / (A * D - B * C);
  r6[py_] = pyf;
  r6[px_] = pxf;
  r6[ct_] = r6[ct_] - (r6[px_] * FX + r6[py_] * FY) * DEL * DEL;
}
