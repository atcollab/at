#include "diff_thinkick.c"

static double B2perp(double bx, double by,
                            double x, double xpr, double y, double ypr)
/* Calculates sqr(|e x B|) , where e is a unit vector in the direction of velocity   */

{
    double v_norm2 = 1.0/(1.0+ SQR(xpr) + SQR(ypr));

	/* components of the velocity vector:
	   ex = xpr;
	   ey = ypr;
	   ez = 1
	*/

	return (SQR(by) + SQR(bx) + SQR(bx*ypr - by*xpr))*v_norm2 ;
}

static void diff_str_kick(double *r6, double *A, double *B, int max_order,
                          double L, double rad_const, double diff_const, double *bdiff)
{
  /* clang-format off */
/*****************************************************************************

The kick is given by

           e L
theta  = - --- B
     x     p    y
            0

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
  double x, xpr, y, ypr, p_norm, B2P, factor;
  double ReSumTemp;

  /* recursively calculate the local transverse magnetic field */
  double ImSum = A[max_order];
  double ReSum = B[max_order];
  for (int i = max_order - 1; i >= 0; i--) {
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

  B2P = B2perp(ImSum, ReSum, x, xpr, y, ypr);
  factor = (1.0  + (SQR(xpr) + SQR(ypr)) / 2.0) / SQR(p_norm) * L;

  if (bdiff) {
    thinkickM(r6, A, B, max_order, L, bdiff);
    thinkickB(r6, ReSum, ImSum, diff_const, B2P, factor, bdiff);
  }

  /* Momentum loss */
  r6[4] -= rad_const * B2P * factor;

  /* recalculate momenta from angles after losing energy */
  p_norm = 1.0 / (1.0 + r6[4]);
  r6[1] = xpr / p_norm;
  r6[3] = ypr / p_norm;

  /* multipole kick */
  r6[1] -= L * ReSum;
  r6[3] += L * ImSum;
}