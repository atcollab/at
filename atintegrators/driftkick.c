/***********************************************************************
 Note: in the US convention the transverse multipole field is written as:

                         max_order+1
                           ----
                           \                       n-1
      (B + iB  )/ B rho  =  >   (ia  + b ) (x + iy)
         y    x            /       n    n
	                       ----
                          n=1
	is a polynomial in (x,y) with the highest order = MaxOrder
	

	Using different index notation 
   
                         max_order
                           ----
                           \                       n
      (B + iB  )/ B rho  =  >   (iA  + B ) (x + iy)
         y    x            /       n    n
	                       ----
                          n=0

	A,B: i=0 ... max_order
   [0] - dipole, [1] - quadrupole, [2] - sextupole ...
   units for A,B[i] = 1/[m]^(i+1)
	Coeficients are stroed in the PolynomA, PolynomB field of the element
	structure in MATLAB

	A[i] (C++,C) =  PolynomA(i+1) (MATLAB) 
	B[i] (C++,C) =  PolynomB(i+1) (MATLAB) 
	i = 0 .. MaxOrder

 ************************************************************************/


static void fastdrift(double* r, double NormL)

/*   NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
 in the loop if momentum deviation (delta) does not change
 such as in 4-th order symplectic integrator w/o radiation
 */

{
   r[0] += NormL*r[1];
   r[2] += NormL*r[3];
   r[5] += NormL*(r[1]*r[1]+r[3]*r[3])/(2*(1+r[4]));
}


static void bndthinkick(double* r, double* A, double* B, double L, double irho, int max_order)
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
   int i;
   double ReSum = B[max_order];
   double ImSum = A[max_order];
   double ReSumTemp;
   /* recursively calculate the local transverse magnetic field
    * Bx = ReSum, By = ImSum
    */
   for (i=max_order-1; i>=0; i--) {
       ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
       ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
       ReSum = ReSumTemp;
   }
   r[1] -=  L*(ReSum-(r[4]-r[0]*irho)*irho);
   r[3] +=  L*ImSum;
   r[5] +=  L*irho*r[0]; /* pathlength */
}


static void strthinkick(double* r, const double* A, const double* B, double L, int max_order)
/***************************************************************************** 
 Calculate and apply a multipole kick to a 6-dimentional
 phase space vector in a straight element (quadrupole)
 
 IMPORTANT !!!
 The reference coordinate system is straight but the field expansion may still
 contain dipole terms: PolynomA(1), PolynomB(1) - in MATLAB notation,
 A[0], B[0] - C,C++ notation
 
 ******************************************************************************/
{
   int i;
   double ReSum = B[max_order];
   double ImSum = A[max_order];
   double ReSumTemp;
   for (i=max_order-1; i>=0; i--) {
      ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
      ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
      ReSum = ReSumTemp;
   }
   r[1] -=  L*ReSum;
   r[3] +=  L*ImSum;
}
