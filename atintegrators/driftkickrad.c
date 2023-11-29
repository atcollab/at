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

static void drift6(double* r, double L)
/*   Input parameter L is the physical length
     1/(1+delta) normalization is done internally
*/
{	double p_norm = 1/(1+r[4]); 
	double NormL  = L*p_norm;   
	r[0]+= NormL*r[1]; 
	r[2]+= NormL*r[3];
	r[5]+= NormL*p_norm*(r[1]*r[1]+r[3]*r[3])/2;
}

#define SQR(X) ((X)*(X))

double StrB2perp(double bx, double by, 
                            double x, double xpr, double y, double ypr)
/* Calculates sqr(|B x e|) , where e is a unit vector in the direction of velocity  */

{	double v_norm2;
	v_norm2 = 1/(1 + SQR(xpr) + SQR(ypr));

	/* components of the normalized velocity vector
	   double ex, ey, ez;
	   ex = xpr; 
	   ey = ypr; 
	   ez = 1;
	*/
  	
	return((SQR(by) + SQR(bx) + SQR(bx*ypr - by*xpr) )*v_norm2) ;

} 


static double B2perp(double bx, double by, double irho,
        double x, double xpr, double y, double ypr)
        /* Calculates sqr(|e x B|) , where e is a unit vector in the direction of velocity  */
        
{
    double v_norm2 = 1/(SQR(1+x*irho)+ SQR(xpr) + SQR(ypr));
    
    /* components of the  velocity vector
     * double ex, ey, ez;
     * ex = xpr;
     * ey = ypr;
     * ez = (1+x*irho);
     */
    
    return((SQR(by*(1+x*irho)) + SQR(bx*(1+x*irho)) + SQR(bx*ypr - by*xpr) )*v_norm2) ;
}

static void bndthinkickrad(double* r, double* A, double* B, double L, double irho, double E0, int max_order)

/*****************************************************************************
Calculate multipole kick in a curved elemrnt (bending magnet)
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
   int i;
   double ImSum = A[max_order];
   double ReSum = B[max_order];
   double ReSumTemp;
   double x ,xpr, y, ypr, p_norm,dp_0, B2P;
   double CRAD = CGAMMA*E0*E0*E0/(TWOPI*1e27);	/* [m]/[GeV^3] M.Sands (4.1) */
   
   /* recursively calculate the local transvrese magnetic field
    * Bx = ReSum, By = ImSum
    */
   for (i=max_order-1; i>=0; i--) {
   	ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
        ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
        ReSum = ReSumTemp;
   }
   
   /* calculate angles from momentums 	*/
   p_norm = 1/(1+r[4]);
   x   = r[0];
   xpr = r[1]*p_norm;
   y   = r[2];
   ypr = r[3]*p_norm;
   
   B2P = B2perp(ImSum, ReSum +irho, irho, x , xpr, y ,ypr);
   
   dp_0 = r[4];
   r[4] = r[4] - CRAD*SQR(1+r[4])*B2P*(1 + x*irho + (SQR(xpr)+SQR(ypr))/2 )*L;
   
   /* recalculate momentums from angles after losing energy for radiation 	*/
   p_norm = 1/(1+r[4]);
   r[1] = xpr/p_norm;
   r[3] = ypr/p_norm;
   
   r[1] -=  L*(ReSum-(dp_0-r[0]*irho)*irho);
   r[3] +=  L*ImSum;
   r[5] +=  L*irho*r[0]; /* pathlength */
}

static void strthinkickrad(double* r, const double* A, const double* B, double L, double E0,int max_order)
/*****************************************************************************
 Calculate and apply a multipole kick to a 6-dimentional
 phase space vector in a straight element ( quadrupole)
 
 IMPORTANT !!!
 he reference coordinate system is straight but the field expansion may still
 ontain dipole terms: PolynomA(1), PolynomB(1) - in MATLAB notation,
 [0], B[0] - C,C++ notation
 
 ******************************************************************************/
{
   int i;
   double ReSum = B[max_order];
   double ImSum = A[max_order];
   double ReSumTemp;
   double irho=0;/*straight elements no curvature.*/
   double x ,xpr, y, ypr, p_norm,dp_0, B2P;
   double CRAD = CGAMMA*E0*E0*E0/(TWOPI*1e27);	/* [m]/[GeV^3] M.Sands (4.1) */
   
   for (i=max_order-1; i>=0; i--) {
      ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
      ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
      ReSum = ReSumTemp;
   }
   
   /* calculate angles from momentums 	*/
   p_norm = 1/(1+r[4]);
   x   = r[0];
   xpr = r[1]*p_norm;
   y   = r[2];
   ypr = r[3]*p_norm;
   
   /*B2P = B2perp(ImSum, ReSum +irho, irho, x , xpr, y ,ypr);*/
   B2P = StrB2perp(ImSum, ReSum , x , xpr, y ,ypr);
   
   dp_0 = r[4];
   r[4] = r[4] - CRAD*SQR(1+r[4])*B2P*(1 + x*irho + (SQR(xpr)+SQR(ypr))/2 )*L;
   
   /* recalculate momentums from angles after losing energy for radiation 	*/
   p_norm = 1/(1+r[4]);
   r[1] = xpr/p_norm;
   r[3] = ypr/p_norm;
   
   r[1] -=  L*(ReSum-(dp_0-r[0]*irho)*irho);
   r[3] +=  L*ImSum;
   r[5] +=  L*irho*r[0]; /* pathlength */
}
