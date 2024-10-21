
static void drift_propagateB(double NormL, double xpr, double ypr, double *bdiff)
{ /* Propagate cumulative Ohmi's diffusion matrix B through a drift.
   B is a (*double) pointer to 1-dimensional array containing 36 elements of
   matrix elements arranged column-by-column as in MATLAB representation

   The relationship between indexes when a 6-by-6 matrix is represented in MATLAB as
   one-dimensional array containing 36 elements arranged column-by-column is
   [i][j] <---> [i+6*j]
  */
  double M66[36];

  /* Initialize M66 to a 6-by-6 identity matrix */
  for (int m = 0; m < 36; m++)
    M66[m] = 0.0;
  for (int m = 0; m < 6; m++)
    M66[m * 7] = 1.0;

  M66[6] = NormL;          /* [0][1] */
  M66[20] = NormL;         /* [2][3] */
  M66[24] = -NormL * xpr;  /* [0][4] */
  M66[26] = -NormL * ypr;  /* [2][4] */
  M66[11] = NormL * xpr;   /* [5][1] */
  M66[23] = NormL * ypr;   /* [5][3] */
  M66[29] = -NormL * (SQR(xpr) + SQR(ypr));   /* [5][4] */

  ATsandwichmmt(M66, bdiff);
}

static void diff_drift(double *r6, double L, double *bdiff)
{
  double p_norm = 1.0 / (1.0+r6[4]);
  double xpr = r6[1] * p_norm;
  double ypr = r6[3] * p_norm;
  double NormL = L * p_norm;

  if (bdiff) {
    drift_propagateB(NormL, xpr, ypr, bdiff);
  }
  r6[0] += NormL * r6[1];
  r6[2] += NormL * r6[3];
  r6[5] += NormL * p_norm * (r6[1]*r6[1] + r6[3]*r6[3]) / 2.0;
}
