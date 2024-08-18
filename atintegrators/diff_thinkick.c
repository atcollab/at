static void thinkickM(const double *r6, const double *A, const double *B, int max_order,
                      double L, double *bdiff)
/* Calculate the symplectic (no radiation) transfer matrix of a
   thin multipole kick near the entrance point r6
*/
{
  double M66[36];
  double ReSumNTemp;
  double ImSumN = max_order * A[max_order];
  double ReSumN = max_order * B[max_order];

  /* Recursively calculate the derivatives
          ReSumN = Re(d(By + iBx)/dx)
          ImSumN = Im(d(By + iBx)/dy)
  */
  for (int n = max_order - 1; n > 0; n--) {
    ReSumNTemp = ReSumN*r6[0] - ImSumN*r6[2] + n*B[n];
    ImSumN = ImSumN*r6[0] + ReSumN*r6[2] + n*A[n];
    ReSumN = ReSumNTemp;
  }

  /* Initialize M66 to a 6-by-6 identity matrix */
  for (int m = 0; m < 36; m++)
    M66[m] = 0.0;
  for (int m = 0; m < 6; m++)
    M66[m*7] = 1.0;

  /*  [i, j] <---> [i+6*j] */
  M66[1] = -L * ReSumN;       /* [1, 0] */
  M66[13] = L * ImSumN;       /* [1, 2] */
  M66[3] = L * ImSumN;        /* [3, 0] */
  M66[15] = L * ReSumN;       /* [3, 2] */

  ATsandwichmmt(M66, bdiff);
}

static void thinkickB(const double *r6, double ReSum, double ImSum,
                      double diff_const, double B2P, double factor, double *bdiff)

/* Calculate Ohmi's diffusion matrix of a thin multipole element */
{
  double B66[36];

  double p_norm = 1.0 / (1.0+r6[4]);
  double p_norm2 = SQR(p_norm);
  double B3P = B2P * sqrt(B2P);
  double BB = diff_const * B3P * factor / p_norm2;  /* m^-1 */

  /* When a 6-by-6 matrix is represented in MATLAB as one-dimentional
     array containing 36 elements arranged column-by-column,
     the relationship between indexes  is
     [i][j] <---> [i+6*j]
  */

  /* initialize B66 to 0 */
  for (int i = 0; i < 36; i++)
    B66[i] = 0.0;

  /* Populate B66 */
  B66[7]  = BB * SQR(r6[1]) * p_norm2;     /* [1][1] */
  B66[19] = BB * r6[1] * r6[3] * p_norm2;  /* [1][3] */
  B66[9]  = B66[19];                       /* [3][1] */
  B66[21] = BB * SQR(r6[3]) * p_norm2;     /* [3][3] */
  B66[10] = BB * r6[1] * p_norm;           /* [4][1] */
  B66[25] = B66[10];                       /* [1][4] */
  B66[22] = BB * r6[3] * p_norm;           /* [4][3] */
  B66[27] = B66[22];                       /* [3][4] */
  B66[28] = BB;                            /* [4][4] */

  ATaddmm(B66, bdiff);
}