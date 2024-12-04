static void thinkickM(const double *r6, double *A, double *B, int max_order,
                      double L, double irho, double *bdiff)
/* Calculate the symplectic (no radiation) transfer matrix of a
   thin multipole kick near the entrance point r6
   For elements with straight coordinate system irho = 0
   For curved elements the B polynomial (PolynomB in MATLAB)
   MUST NOT include the guide field  By0 = irho * E0 /(c*e)
*/
{
  double M66[36];
  double ReSumNTemp;
  double ImSumN = max_order * A[max_order];
  double ReSumN = max_order * B[max_order];

  /* Recursively calculate the derivatives
          ReSumN = (irho/B0)*Re(d(By + iBx)/dx)
          ImSumN = (irho/B0)*Im(d(By + iBx)/dy)
   */
  for (int n = max_order - 1; n > 0; n--) {
    ReSumNTemp = (ReSumN * r6[0] - ImSumN * r6[2]) + n * B[n];
    ImSumN = ImSumN * r6[0] + ReSumN * r6[2] + n * A[n];
    ReSumN = ReSumNTemp;
  }

  /* Initialize M66 to a 6-by-6 identity matrix */
  for (int m = 0; m < 36; m++)
    M66[m] = 0.0;
  for (int m = 0; m < 6; m++)
    M66[m*7] = 1.0;

  /* The relationship between indexes when a 6-by-6 matrix is
     represented in MATLAB as one-dimentional array containing
     36 elements arranged column-by-column is
     [i][j] <---> [i+6*j]
  */

  M66[1] = -L * ReSumN;       /* [1][0] */
  M66[13] = L * ImSumN;       /* [1][2] */
  M66[3] = L * ImSumN;        /* [3][0] */
  M66[15] = L * ReSumN;       /* [3][2] */
  M66[25] = L * irho;         /* [1][4] */
  M66[1] += -L * irho * irho; /* [1][0] */
  M66[5] = L * irho;          /* [5][0] */

  ATsandwichmmt(M66, bdiff);
}

static void thinkickB(double *r6, double ReSum, double ImSum,
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

static double B2perp(double bx, double by, double irho,
                            double x, double xpr, double y, double ypr)
/* Calculates sqr(|e x B|) , where e is a unit vector in the direction of velocity   */

{
    double v_norm2 = 1.0/(SQR(1.0+x*irho)+ SQR(xpr) + SQR(ypr));

	/* components of the  velocity vector:
	   ex = xpr;
	   ey = ypr;
	   ez = (1+x*irho);
	*/

	return (SQR(by*(1+x*irho)) + SQR(bx*(1+x*irho)) + SQR(bx*ypr - by*xpr))*v_norm2 ;
}

static void diff_bnd_kick(double *r6, double *A, double *B, int max_order,
                          double L, double irho, double rad_const, double diff_const, double *bdiff) {
  /* clang-format off */
/*****************************************************************************
The design magnetic field Byo that provides this curvature By0 = irho * E0 /(c*e)
MUST NOT be included in the dipole term PolynomB(1)(MATLAB notation)(B[0] C notation)
of the By field expansion
HOWEVER!!! to calculate the effect of classical radiation the full field must be
used in the square of the |v x B|.
When calling B2perp(Bx, By, ...), use the By = ReSum + irho, where ReSum is the
normalized vertical field - sum of the polynomial terms in PolynomB.

The kick is given by

           e L      L delta      L x
theta  = - --- B  + -------  -  -----  ,
     x     p    y     rho           2
            0                    rho

         e L
theta  = --- B
     y    p   x
           0

                          max_order
                            ----
                            \                       n
	   (B + iB  )/ B rho  =  >   (iA  + B ) (x + iy)
         y    x             /       n    n
  	                        ----
                            n=0

  ******************************************************************************/
  /* clang-format on */
  int i;
  double ImSum = A[max_order];
  double ReSum = B[max_order];
  double x, xpr, y, ypr, p_norm, dp_0, B2P, factor;
  double ReSumTemp;

  /* recursively calculate the local transverse magnetic field */
  for (i = max_order - 1; i >= 0; i--) {
    ReSumTemp = ReSum * r6[0] - ImSum * r6[2] + B[i];
    ImSum = ImSum * r6[0] + ReSum * r6[2] + A[i];
    ReSum = ReSumTemp;
  }
  /* calculate angles from momenta */
  p_norm = 1.0 / (1.0+r6[4]);
  x = r6[0];
  xpr = r6[1] * p_norm;
  y = r6[2];
  ypr = r6[3] * p_norm;

  B2P = B2perp(ImSum, ReSum+irho, irho, x, xpr, y, ypr);
  factor = (1.0 + x*irho + (SQR(xpr) + SQR(ypr)) / 2.0) / SQR(p_norm) * L;

  if (bdiff) {
    thinkickM(r6, A, B, max_order, L, irho, bdiff);
    thinkickB(r6, ReSum, ImSum, diff_const, B2P, factor, bdiff);
  }

  dp_0 = r6[4]; /* save a copy of the initial value of dp/p */

  r6[4] -= rad_const * B2P * factor;

  /* recalculate momenta from angles after losing energy */
  p_norm = 1.0 / (1.0 + r6[4]);
  r6[1] = xpr / p_norm;
  r6[3] = ypr / p_norm;

  r6[1] -= L * (ReSum - (dp_0 - r6[0] * irho) * irho);
  r6[3] += L * ImSum;
  r6[5] += L * irho * r6[0]; /* pathlength */
}