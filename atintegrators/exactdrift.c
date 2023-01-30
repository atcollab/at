
#define SQR(X) ((X)*(X))

static double get_pz(double *r_in) {
  return sqrt(SQR(1 + r_in[4]) - SQR(r_in[1]) - SQR(r_in[3]));
}

static void exact_drift(double *r_in, double L) {
  double NormL = L / get_pz(r_in);
  r_in[0] += r_in[1] * NormL;
  r_in[2] += r_in[3] * NormL;
  r_in[5] += NormL * (1.0 + r_in[4]) - L;
}
