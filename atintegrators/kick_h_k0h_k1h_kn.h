/***********************************************************************
Expansion of the magnetic field in AT:

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
*************************************************************************/

#ifdef RADIATION
#include "diff_thinkick.c"

static double B2perp(double bx, double by, double irho, double x, double xpr, double y, double ypr)
/* Calculates sqr(|B x e|), where e is a unit vector in the direction of velocity */
{
    double xh1 = SQR(1.0 + x*irho);
    double v_norm2 = (xh1 + SQR(xpr) + SQR(ypr));

	/* components of the velocity vector:
	   ex = xpr;
	   ey = ypr;
	   ez = (1+x*irho);

	   sqr(|B x e|) = sqr(|B|) * sqr(|e|) - sqr(B.e)
	*/

    return((xh1*SQR(by) + xh1*SQR(bx) + SQR(bx*ypr - by*xpr)) / v_norm2) ;
}

static void kick(double *r6, double A0, double B0, const double *A, const double *B, int max_order,
                 double L, double irho, double rad_const, double diff_const, double *bdiff)
#else
static void kick(double *r6, double A0, double B0, const double *A, const double *B, int max_order,
                 double L, double irho)
#endif /* RADIATION */
{
/* clang-format off */
/*****************************************************************************
Calculate multipole kick in a curved element (bending magnet)
The reference coordinate system has the curvature given by the inverse (design) radius irho.
IMPORTANT !!!
The magnetic field Bo that provides this curvature MUST NOT be included in the dipole term

The kick is given by
                                                 2   2
           e L      L delta      L x      L K1 (x - y /2)
theta  = - --- B  + -------  -  -----  - -----------------,
     x     p    y     rho           2           rho
            0                    rho

         e L       L K1 x y
theta  = --- B  + ----------
     y    p   x      rho
           0

******************************************************************************/
/* clang-format on */
    double ReSum = B[max_order];
    double ImSum = A[max_order];
    double ReSumTemp;
    double x = r6[0];
    double y = r6[2];
    double dp_0 = r6[4]; /* save a copy of the initial value of dp/p */
    double B1 = (max_order >= 1) ? B[1] : 0.0;

    /* recursively calculate the local transverse magnetic field */
    for (int i = max_order - 1; i >= 0; i--) {
        ReSumTemp = ReSum*x - ImSum*y + B[i];
        ImSum = ImSum*x + ReSum*y + A[i];
        ReSum = ReSumTemp;
    }
    ReSum += B0;
    ImSum += A0;

    #ifdef RADIATION
    double p_norm = 1.0 / (1.0+r6[4]);

    /* calculate angles from momenta */
    double xpr = r6[1] * p_norm;
    double ypr = r6[3] * p_norm;

    #ifdef CURVATURE_IN_B0
    double B2P = B2perp(ImSum, ReSum, irho, x, xpr, y ,ypr);
    #else
    double B2P = B2perp(ImSum, ReSum + irho, irho, x, xpr, y ,ypr);
    #endif
    double factor = L * (1.0 + x*irho + (SQR(xpr) + SQR(ypr)) / 2.0) / SQR(p_norm);

    /* Propagation of the diffusion matrix */
    if (bdiff) {
        thinkickM(r6, A, B, max_order, L, irho, bdiff);
        thinkickB(r6, ReSum, ImSum, diff_const, B2P, factor, bdiff);
    }

    /* Momentum loss */
    r6[4] -= rad_const * B2P * factor;

    /* Recalculate momenta from angles after losing energy */
    p_norm = 1.0 / (1.0 + r6[4]);
    r6[1] = xpr / p_norm;
    r6[3] = ypr / p_norm;
    #endif /* RADIATION */

    /* Multipole kick */
    r6[1] -= L * (ReSum + irho*(x*irho - dp_0 + B1*(x*x-0.5*y*y)));
    r6[3] += L * (ImSum + irho*B1*x*y);
    r6[5] += L * irho*x; /* pathlength */
}
