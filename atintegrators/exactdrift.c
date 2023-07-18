#include <math.h>
#define SQR(X) ((X)*(X))

static double get_pz(double *r6) {
  return sqrt(SQR(1 + r6[4]) - SQR(r6[1]) - SQR(r6[3]));
}

/* Forest 10.23, exact drift
   L: length [m]
*/
static void exact_drift(double *r6, double L) {
  double NormL = L / get_pz(r6);
  r6[0] += r6[1] * NormL;
  r6[2] += r6[3] * NormL;
  r6[5] += NormL * (1.0 + r6[4]);   /* Absolute path length */
}
