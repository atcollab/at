/***********************************************************************
 Note: in the US convention the transverse multipole field is written as:

                      max_order+1
                       ---
                       \                       n-1
  (B + iB  )/ B rho  =  >   (ia  + b ) (x + iy)
    y    x             /       n    n
                       ----
                       n=1
 is a polynomial in (x,y) with the highest order = MaxOrder


 Using different index notation

                      max_order
                      ----
                      \                       n
 (B + iB  )/ B rho  =  >   (iA  + B ) (x + iy)
   y    x             /       n    n
                      ----
                      n=0

 A,B: i=0 ... max_order
 [0] - dipole, [1] - quadrupole, [2] - sextupole ...
 units for A,B[i] = 1/[m]^(i+1)
 Coefficients are stored in the PolynomA, PolynomB field of the element
 structure in MATLAB

 A[i] (C++,C) =  PolynomA(i+1) (MATLAB)
 B[i] (C++,C) =  PolynomB(i+1) (MATLAB)
 i = 0 .. MaxOrder

*************************************************************************/

#define SQR(X) ((X)*(X))

static double B2perp(double bx, double by, double irho, double x, double xpr, double y, double ypr)
/* Calculates sqr(|e x B|), where e is a unit vector in the direction of velocity */
{
    /* components of the  velocity vector
       double ex, ey, ez;
       ex = xpr;
       ey = ypr;
       ez = (1+x*irho) * sqrt(1 - xpr^2 - ypr^2);

	   sqr(|B x e|) = sqr(|B|) * sqr(|e|) - sqr(B.e)
    */
    double nrm = SQR(1.0+x*irho);
    double v_norm2 = nrm + SQR(xpr)*(1.0-nrm) + SQR(ypr)*(1.0-nrm);

    return SQR(bx) + SQR(by) - SQR(bx*xpr + by*ypr)/v_norm2;
}

//static void ex_bndthinkickrad(double* r, double* A, double* B, double L, double irho, double E0, int max_order)
static void kick(double* r6, double A0, double B0, double* A, double* B, int max_order,
                              double L, double irho, double rad_const, double diff_const, double *bdiff)

/*****************************************************************************
Calculate multipole kick in a curved element (bending magnet)
The reference coordinate system  has the curvature given by the inverse
(design) radius irho.
IMPORTANT !!!
The magnetic field Bo that provides this curvature MUST NOT be included in the dipole term
PolynomB[1](MATLAB notation)(C: B[0] in this function) of the By field expansion
HOWEVER!!! to calculate the effect of classical radiation the full field must be
used in the square of the |v x B|.
When calling B2perp(Bx, By, ...), use the By = RESum + irho, where ImSum is the sum of
the polynomial terms in PolynomB.

 The kick is given by

             e L      L delta      L x
  theta  = - --- B  + -------  -  -----  ,
       x      p   y     rho           2
               0                   rho

           e L
  theta  = --- B
       y    p   x
             0

 ******************************************************************************/
{
  double ReSum = B[max_order];
  double ImSum = A[max_order];
  double ReSumTemp;
  double B2P;
  double p_norm = 1.0 / (1.0+r6[4]);
  double x = r6[0];
  double y = r6[2];
  double dp_0 = r6[4]; /* save a copy of the initial value of dp/p */

  /* calculate angles from momenta */
  double xpr = r6[1] * p_norm;
  double ypr = r6[3] * p_norm;

  /* recursively calculate the local transverse magnetic field */
  for (int i = max_order - 1; i >= 0; i--) {
    ReSumTemp = ReSum*x - ImSum*y + B[i];
    ImSum = ImSum*x + ReSum*y + A[i];
    ReSum = ReSumTemp;
  }
  ReSum += B0;
  ImSum += A0;

  B2P = B2perp(ImSum, ReSum+irho, irho, x , xpr, y ,ypr);

  /* Momentum loss */
  r6[4] -= rad_const * SQR(1+r6[4]) * B2P * (1.0+x*irho) * L / sqrt(1.0 - xpr*xpr - ypr*ypr);

  /* recalculate momentums from angles after losing energy for radiation 	*/
  p_norm = 1.0 / (1.0+r6[4]);
  r6[1] = xpr/p_norm;
  r6[3] = ypr/p_norm;

  /* Multipole kick */
  r6[1] -= L * (ReSum + irho*B[1]*x*y);
  r6[3] += L * (ImSum - irho*B[1]*(x*x-y*y/2.0));
}
