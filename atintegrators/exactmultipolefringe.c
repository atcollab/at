#define x_ 0
#define px_ 1
#define y_ 2
#define py_ 3
#define delta_ 4
#define ct_ 5

static void multipole_fringe(double *r_in, double L,
                             double *polya, double *polyb, int max_order,
                             int edge)
{
  // PTC multipole_fringer
  // Forest 13.29
  // not re-derived and checked
  // note this is the sum over n of Forest 13.29
  // one for each multipole component

  double I, U, V, DU, DV, DUX, DVX, DUY, DVY, FX, FY, FX_X, FX_Y, FY_X, FY_Y,
      RX, IX, DRX, DIX;

  if (edge == 0) {
    I = 1;
  } else {
    I = -1;
  }

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

    RX = DRX * r_in[x_] - DIX * r_in[y_];
    IX = DRX * r_in[y_] + DIX * r_in[x_];

    U = B * RX - A * IX;
    V = B * IX + A * RX;
    DU = B * DRX - A * DIX;
    DV = B * DIX + A * DRX;

    double f1 = -I / 4.0 / (j + 1);

    U = U * f1;
    V = V * f1;
    DU = DU * f1;
    DV = DV * f1;

    DUX = j * DU;
    DVX = j * DV;
    DUY = -j * DV;
    DVY = j * DU;

    double nf = 1.0 * (j + 2) / j;

    FX += U * r_in[x_] + nf * V * r_in[y_];
    FY += U * r_in[y_] - nf * V * r_in[x_];

    FX_X += DUX * r_in[x_] + U + nf * r_in[y_] * DVX;
    FX_Y += DUY * r_in[x_] + nf * V + nf * r_in[y_] * DVY;

    FY_X += DUX * r_in[y_] - nf * V - nf * r_in[x_] * DVX;
    FY_Y += DUY * r_in[y_] + U - nf * r_in[x_] * DVY;
  }

  double DEL = 1.0 / (1 + r_in[delta_]);

  // solve 2x2 matrix equation

  double A = 1 - FX_X * DEL;
  double B = -FY_X * DEL;
  double D = 1 - FY_Y * DEL;
  double C = -FX_Y * DEL;

  r_in[x_] = r_in[x_] - FX * DEL;
  r_in[y_] = r_in[y_] - FY * DEL;

  double pxf = (D * r_in[px_] - B * r_in[py_]) / (A * D - B * C);
  double pyf = (A * r_in[py_] - C * r_in[px_]) / (A * D - B * C);
  r_in[py_] = pyf;
  r_in[px_] = pxf;
  r_in[ct_] = r_in[ct_] - (r_in[px_] * FX + r_in[py_] * FY) * DEL * DEL;
}
