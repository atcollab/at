#include "diff_thinkick.c"

static double B2perp(double bx, double by,
                            double x, double xpr, double y, double ypr)
/* Calculates sqr(|B x e|) , where e is a unit vector in the direction of velocity  */

{
    /* components of the normalized velocity vector
	   double ex, ey, ez;
	   ex = xpr;
	   ey = ypr;
	   ez = sqrt(1 - xpr^2 - ypr^2);

	   sqr(|B x e|) = sqr(|B|) * sqr(|e|) - sqr(B.e)
	*/

  	return SQR(bx) + SQR(by) - SQR(bx*xpr + by*ypr);
}

static void diff_str_exactkick(double* r6, const double* A, const double* B, int max_order,
    double L, double rad_const, double diff_const, double *bdiff)
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
{
  double x ,xpr, y, ypr, p_norm, B2P, factor;
  double ReSumTemp;

  /* recursively calculate the local transverse magnetic field */
  double ReSum = B[max_order];
  double ImSum = A[max_order];
  for (int i=max_order-1; i>=0; i--) {
    ReSumTemp = ReSum*r6[0] - ImSum*r6[2] + B[i];
    ImSum = ImSum*r6[0] + ReSum*r6[2] + A[i];
    ReSum = ReSumTemp;
  }

  /* calculate angles from momentums 	*/
  p_norm = 1.0 / (1.0+r6[4]);
  x = r6[0];
  xpr = r6[1] * p_norm;
  y = r6[2];
  ypr = r6[3] * p_norm;

  B2P = B2perp(ImSum, ReSum, x, xpr, y ,ypr);
  factor = L / SQR(p_norm) / sqrt(1.0 - SQR(xpr) - SQR(ypr));

  if (bdiff) {
    thinkickM(r6, A, B, max_order, L, bdiff);
    thinkickB(r6, ReSum, ImSum, diff_const, B2P, factor, bdiff);
  }

  /* Momentum loss */
  r6[4] -= rad_const * B2P * factor;

  /* recalculate momentums from angles after losing energy for radiation 	*/
  p_norm = 1.0 / (1.0+r6[4]);
  r6[1] = xpr / p_norm;
  r6[3] = ypr / p_norm;

  /* multipole kick */
  r6[1] -= L * ReSum;
  r6[3] += L * ImSum;
}
