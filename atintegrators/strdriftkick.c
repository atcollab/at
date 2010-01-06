
static void strthinkick(double* r, const double* A, const double* B, double L, int max_order)
/***************************************************************************** 
 Calculate and apply a multipole kick to a 6-dimentional
 phase space vector in a straight element ( quadrupole)
 
 IMPORTANT !!!
 The reference coordinate system is straight but the field expansion may still
 contain dipole terms: PolynomA(1), PolynomB(1) - in MATLAB notation,
 A[0], B[0] - C,C++ notation
 
 
 Note: in the US convention the transverse multipole field is written as:
 
                      max_order+1
                       ---
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
 
 ******************************************************************************/
{  int i;
   double ReSum = B[max_order];
   double ImSum = A[max_order];
   double ReSumTemp;
   for(i=max_order-1;i>=0;i--)
   {   ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
      ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
      ReSum = ReSumTemp;
   }
   
   r[1] -=  L*ReSum;
   r[3] +=  L*ImSum;
}

static void fastdrift(double* r, double NormL)

/*   NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
 in the loop if momentum deviation (delta) does not change
 such as in 4-th order symplectic integrator w/o radiation
 */

{   double dx = NormL*r[1];
   double dy = NormL*r[3];
   r[0]+= dx;
   r[2]+= dy;
   r[5]+= NormL*(r[1]*r[1]+r[3]*r[3])/(2*(1+r[4]));
}
