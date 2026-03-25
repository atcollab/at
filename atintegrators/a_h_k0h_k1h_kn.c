static void kick(double* r, double A0, double B0, double* A, double* B, double L, double irho, int max_order)
/*****************************************************************************
Calculate multipole kick in a curved elemrnt (bending magnet)
The reference coordinate system  has the curvature given by the inverse
(design) radius irho.
IMPORTANT !!!
The magnetic field Bo that provides this curvature MUST NOT be included in the dipole term
PolynomB[1](MATLAB notation)(C: B[0] in this function) of the By field expansion

The kick is given by

           e L      L delta      L x
theta  = - --- B  + -------  -  -----  ,
     x     p    y     rho           2
            0                    rho

         e L
theta  = --- B
     y    p   x
           0

*************************************************************************/
{
   double ReSum = B[max_order];
   double ImSum = A[max_order];
   double ReSumTemp;
   double x = r[0];
   double y = r[2];

   /* recursively calculate the local transverse magnetic field */
   for (int i=max_order-1; i>=0; i--) {
       ReSumTemp = ReSum*x - ImSum*y + B[i];
       ImSum = ImSum*x +  ReSum*y + A[i];
       ReSum = ReSumTemp;
   }
   ReSum += B0;
   ImSum += A0;

   r[1] -=  L*(ReSum - (r[4]-x*irho)*irho + irho*B[1]*x*y);
   r[3] +=  L*(ImSum - irho*B[1]*(x*x-y*y/2.0));
   r[5] +=  L*irho*x; /* pathlength */
}
