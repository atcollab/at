#ifndef PXYZ
#define PXYZ
static double pxyz(double dp1, double px, double py)
{
  return sqrt(dp1*dp1 - px*px - py*py);
}
#endif /*PXYZ*/

static yrot_propagate(double *r6, double c, double s, double pz, double p, double *bdiff)
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

    if (phi != 0.0) {
        double dp1 = 1.0 + r6[delta_];
        double c = cos(phi);
        double s = sin(phi);
        double pz = pxyz(dp1, r6[px_], r6[py_]);
        double p = c*pz - s*r6[px_];
        double px = s*pz + c*r6[px_];
        double x = r6[x_]*pz/p;
        double dy = r6[x_]*r6[py_]*s/p;
        double dct = dp1*r6[x_]*s/p;

        if (bdiff) {
            yrot_propagate(r6, c, s, pz, p, bdiff);
        }
        r6[x_] = x;
        r6[px_] = px;
        r6[y_] += dy;
        r6[ct_] += dct;
    }
}