#include <math.h>
#include "atlalib.c"
/*
To speed up the integration loop, the path length is computed in absolute on
each step, and the reference total length is subtracted at the end of the loop.
For other uses, the relative path length is computed.
*/
#ifdef MAGNET_PASS
#define ABSOLUTE_PATH_LENGTH
#define FIX_LENGTH(length) r6[5] -= (length)
#endif

#define SQR(X) ((X)*(X))

#define DRIFT(r6, length, irho, bdiff) drift(r6, length, bdiff)

static void drift_propagateB(double NormL, double xpr, double ypr, double *bdiff)
{ /* Propagate cumulative Ohmi's diffusion matrix B through an exact drift.
   B is a (*double) pointer to 1-dimensional array containing 36 elements of
   matrix elements arranged column-by-column as in MATLAB representation

   The relationship between indexes when a 6-by-6 matrix is represented in MATLAB as
   one-dimensional array containing 36 elements arranged column-by-column is
   [i][j] <---> [i+6*j]
  */
  double M66[36];
  double den = 1.0 - SQR(xpr) - SQR(ypr);

  /* Initialize M66 to a 6-by-6 identity matrix */
  for (int m = 0; m < 36; m++)
    M66[m] = 0.0;
  for (int m = 0; m < 6; m++)
    M66[m * 7] = 1.0;

  M66[6]  = NormL * (1.0 - SQR(ypr)) / den;    /* [0, 1] */
  M66[18] = NormL * xpr * ypr / den;           /* [0, 3] */
  M66[24] = -NormL * xpr / den;                /* [0, 4] */
  M66[8]  = NormL * xpr * ypr / den;           /* [2, 1] */
  M66[20] = NormL * (1.0 - SQR(xpr)) / den;    /* [2, 3] */
  M66[26] = -NormL * ypr / den;                /* [2, 4] */
  M66[11] = NormL * xpr / den;                 /* [5, 1] */
  M66[23] = NormL * ypr / den;                 /* [5, 3] */
  M66[29] = -NormL * (SQR(xpr) + SQR(ypr)) / den;   /* [5][4] */

  ATsandwichmmt(M66, bdiff);
}
/* Forest 10.23, exact drift
   L: length [m]
*/
static void drift(double *r6, double L, double *bdiff)
{
  double p_norm = 1.0 / (1.0+r6[4]);
  double xpr = r6[1] * p_norm;
  double ypr = r6[3] * p_norm;
  double NormL = L * p_norm / sqrt(1.0 - SQR(xpr) - SQR(ypr));

  if (bdiff) {
    drift_propagateB(NormL, xpr, ypr, bdiff);
  }

  r6[0] += NormL * r6[1];
  r6[2] += NormL * r6[3];
  r6[5] += NormL * (1.0 + r6[4]);   /* Absolute path length */
  #ifndef ABSOLUTE_PATH_LENGTH
  r6[5] -= L;
  #endif
}
