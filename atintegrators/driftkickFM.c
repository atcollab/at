#define TWOPI     6.28318530717959
#define CGAMMA 	  8.846056192e-05 			/* [m]/[GeV^3] Ref[1] (4.1)      */

#define SQR(X) ((X)*(X))

struct elem
{
    double Length;
    int NumIntSteps;
    double *LUT_Bx;
    double *LUT_By;
    double *LUT_dBxdx;
    double *LUT_dBydx;
    double *LUT_dBxdy;
    double *LUT_dBydy;
    double *LUT_d2Bxdxdx;
    double *LUT_d2Bxdxdy;
    double *LUT_d2Bxdydy;
    double *LUT_d2Bydxdx;
    double *LUT_d2Bydxdy;
    double *LUT_d2Bydydy;
    double *X;
    double *Y;
    int Nx;
    int Ny;
    double xmin;
    double ymin;
    double delta_x;
    double delta_y;
    double Energy;
    /* Optional fields */
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *KickAngle;
};


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

void get_xy_indeces1(double x,double y,int Nx,int Ny,double xmin,double ymin,double delta_x,double delta_y,int *xi,int *yi,double *dx,double *dy)
{
    /*printf("here I'm finding the indeces for \n");
    printf("x = %f and y = %f \n",x,y);*/
    /**xi=(int)floor(0.5+((x-xmin)/delta_x));
    *yi=(int)floor(0.5+((y-ymin)/delta_y));*/
    *xi=(int)floor((x-xmin)/delta_x);
    *yi=(int)floor((y-ymin)/delta_y);
    if(*xi<1) *xi=1;
    if(*xi>Nx-2) *xi=Nx-2;
    if(*yi<1) *yi=1;
    if(*yi>Ny-2) *yi=Ny-2;
    *dx = x-*xi*delta_x-xmin;
    *dy = y-*yi*delta_y-ymin;
}
void get_xy_indeces2(double x,double y,int Nx,int Ny,double xmin,double ymin,double delta_x,double delta_y,int *xi,int *yi,double *dx,double *dy)
{
    /*printf("here I'm finding the indeces for \n");
    printf("x = %f and y = %f \n",x,y);*/
    /**xi=(int)floor(0.5+((x-xmin)/delta_x));
    *yi=(int)floor(0.5+((y-ymin)/delta_y));*/
    
    *xi=(int)floor((x-xmin)/delta_x);
    *yi=(int)floor((y-ymin)/delta_y);
    if(*xi<1) *xi=1;
    if(*xi>Nx-2) *xi=Nx-2;
    if(*yi<1) *yi=1;
    if(*yi>Ny-2) *yi=Ny-2;
    *dx = x-*xi*delta_x-xmin;
    *dy = y-*yi*delta_y-ymin;
}

void interpolate2ord(struct elem *Elem,int xi,int yi,double dx,double dy,double *Bx,double *By)
{
    int Nx=Elem->Nx;
    int Ny=Elem->Ny;
    
    double Bx_sw, By_sw, Bx_se, By_se, Bx_nw, By_nw, Bx_ne, By_ne;
    
    double Bx0_sw, dBxdx_sw, dBxdy_sw, d2Bxdxdx_sw, d2Bxdxdy_sw, d2Bxdydy_sw;
    double By0_sw, dBydx_sw, dBydy_sw, d2Bydxdx_sw, d2Bydxdy_sw, d2Bydydy_sw;
    double Bx0_se, dBxdx_se, dBxdy_se, d2Bxdxdx_se, d2Bxdxdy_se, d2Bxdydy_se;
    double By0_se, dBydx_se, dBydy_se, d2Bydxdx_se, d2Bydxdy_se, d2Bydydy_se;
    double Bx0_nw, dBxdx_nw, dBxdy_nw, d2Bxdxdx_nw, d2Bxdxdy_nw, d2Bxdydy_nw;
    double By0_nw, dBydx_nw, dBydy_nw, d2Bydxdx_nw, d2Bydxdy_nw, d2Bydydy_nw;
    double Bx0_ne, dBxdx_ne, dBxdy_ne, d2Bxdxdx_ne, d2Bxdxdy_ne, d2Bxdydy_ne;
    double By0_ne, dBydx_ne, dBydy_ne, d2Bydxdx_ne, d2Bydxdy_ne, d2Bydydy_ne;
    
    double dx_sw=dx, dy_sw=dy, dx2_sw=dx*dx, dxdy_sw=dx*dy, dy2_sw=dy*dy;
    double dx_se=dx-Elem->delta_x, dy_se=dy, dx2_se=(dx-Elem->delta_x)*(dx-Elem->delta_x), dxdy_se=(dx-Elem->delta_x)*dy, dy2_se=dy*dy;
    double dx_nw=dx, dy_nw=dy-Elem->delta_y, dx2_nw=dx*dx, dxdy_nw=dx*(dy-Elem->delta_y), dy2_nw=(dy-Elem->delta_y)*(dy-Elem->delta_y);
    double dx_ne=dx-Elem->delta_x, dy_ne=dy-Elem->delta_y, dx2_ne=(dx-Elem->delta_x)*(dx-Elem->delta_x), dxdy_ne=(dx-Elem->delta_x)*(dy-Elem->delta_y), dy2_ne=(dy-Elem->delta_y)*(dy-Elem->delta_y);
    
    double deltaxdeltay=(Elem->delta_x) * (Elem->delta_y);

    /*dx2=dx*dx;
    dxdy=dx*dy;
    dy2=dy*dy;*/
    
    Bx0_sw=Elem->LUT_Bx[yi*Nx+xi];
    dBxdx_sw=Elem->LUT_dBxdx[yi*Nx+xi];
    dBxdy_sw=Elem->LUT_dBxdy[yi*Nx+xi];
    d2Bxdxdx_sw=Elem->LUT_d2Bxdxdx[yi*Nx+xi];
    d2Bxdxdy_sw=Elem->LUT_d2Bxdxdy[yi*Nx+xi];
    d2Bxdydy_sw=Elem->LUT_d2Bxdydy[yi*Nx+xi];
    By0_sw=Elem->LUT_By[yi*Nx+xi];
    dBydx_sw=Elem->LUT_dBydx[yi*Nx+xi];
    dBydy_sw=Elem->LUT_dBydy[yi*Nx+xi];
    d2Bydxdx_sw=Elem->LUT_d2Bydxdx[yi*Nx+xi];
    d2Bydxdy_sw=Elem->LUT_d2Bydxdy[yi*Nx+xi];
    d2Bydydy_sw=Elem->LUT_d2Bydydy[yi*Nx+xi];
    
    Bx0_se=Elem->LUT_Bx[yi*Nx+xi+1];
    dBxdx_se=Elem->LUT_dBxdx[yi*Nx+xi+1];
    dBxdy_se=Elem->LUT_dBxdy[yi*Nx+xi+1];
    d2Bxdxdx_se=Elem->LUT_d2Bxdxdx[yi*Nx+xi+1];
    d2Bxdxdy_se=Elem->LUT_d2Bxdxdy[yi*Nx+xi+1];
    d2Bxdydy_se=Elem->LUT_d2Bxdydy[yi*Nx+xi+1];
    By0_se=Elem->LUT_By[yi*Nx+xi+1];
    dBydx_se=Elem->LUT_dBydx[yi*Nx+xi+1];
    dBydy_se=Elem->LUT_dBydy[yi*Nx+xi+1];
    d2Bydxdx_se=Elem->LUT_d2Bydxdx[yi*Nx+xi+1];
    d2Bydxdy_se=Elem->LUT_d2Bydxdy[yi*Nx+xi+1];
    d2Bydydy_se=Elem->LUT_d2Bydydy[yi*Nx+xi+1];
    
    Bx0_nw=Elem->LUT_Bx[(yi+1)*Nx+xi];
    dBxdx_nw=Elem->LUT_dBxdx[(yi+1)*Nx+xi];
    dBxdy_nw=Elem->LUT_dBxdy[(yi+1)*Nx+xi];
    d2Bxdxdx_nw=Elem->LUT_d2Bxdxdx[(yi+1)*Nx+xi];
    d2Bxdxdy_nw=Elem->LUT_d2Bxdxdy[(yi+1)*Nx+xi];
    d2Bxdydy_nw=Elem->LUT_d2Bxdydy[(yi+1)*Nx+xi];
    By0_nw=Elem->LUT_By[(yi+1)*Nx+xi];
    dBydx_nw=Elem->LUT_dBydx[(yi+1)*Nx+xi];
    dBydy_nw=Elem->LUT_dBydy[(yi+1)*Nx+xi];
    d2Bydxdx_nw=Elem->LUT_d2Bydxdx[(yi+1)*Nx+xi];
    d2Bydxdy_nw=Elem->LUT_d2Bydxdy[(yi+1)*Nx+xi];
    d2Bydydy_nw=Elem->LUT_d2Bydydy[(yi+1)*Nx+xi];
    
    
    Bx0_ne=Elem->LUT_Bx[(yi+1)*Nx+xi+1];
    dBxdx_ne=Elem->LUT_dBxdx[(yi+1)*Nx+xi+1];
    dBxdy_ne=Elem->LUT_dBxdy[(yi+1)*Nx+xi+1];
    d2Bxdxdx_ne=Elem->LUT_d2Bxdxdx[(yi+1)*Nx+xi+1];
    d2Bxdxdy_ne=Elem->LUT_d2Bxdxdy[(yi+1)*Nx+xi+1];
    d2Bxdydy_ne=Elem->LUT_d2Bxdydy[(yi+1)*Nx+xi+1];
    By0_ne=Elem->LUT_By[(yi+1)*Nx+xi+1];
    dBydx_ne=Elem->LUT_dBydx[(yi+1)*Nx+xi+1];
    dBydy_ne=Elem->LUT_dBydy[(yi+1)*Nx+xi+1];
    d2Bydxdx_ne=Elem->LUT_d2Bydxdx[(yi+1)*Nx+xi+1];
    d2Bydxdy_ne=Elem->LUT_d2Bydxdy[(yi+1)*Nx+xi+1];
    d2Bydydy_ne=Elem->LUT_d2Bydydy[(yi+1)*Nx+xi+1];

    Bx_sw=Bx0_sw+dBxdx_sw*dx_sw+dBxdy_sw*dy_sw+(d2Bxdxdx_sw*dx2_sw+d2Bxdxdy_sw*dxdy_sw+d2Bxdydy_sw*dy2_sw);
    By_sw=By0_sw+dBydx_sw*dx_sw+dBydy_sw*dy_sw+(d2Bydxdx_sw*dx2_sw+d2Bydxdy_sw*dxdy_sw+d2Bydydy_sw*dy2_sw);
    
    Bx_se=Bx0_se+dBxdx_se*dx_se+dBxdy_se*dy_se+(d2Bxdxdx_se*dx2_se+d2Bxdxdy_se*dxdy_se+d2Bxdydy_se*dy2_se);
    By_se=By0_se+dBydx_se*dx_se+dBydy_se*dy_se+(d2Bydxdx_se*dx2_se+d2Bydxdy_se*dxdy_se+d2Bydydy_se*dy2_se);
    
    Bx_nw=Bx0_nw+dBxdx_nw*dx_nw+dBxdy_nw*dy_nw+(d2Bxdxdx_nw*dx2_nw+d2Bxdxdy_nw*dxdy_nw+d2Bxdydy_nw*dy2_nw);
    By_nw=By0_nw+dBydx_nw*dx_nw+dBydy_nw*dy_nw+(d2Bydxdx_nw*dx2_nw+d2Bydxdy_nw*dxdy_nw+d2Bydydy_nw*dy2_nw);
    
    Bx_ne=Bx0_ne+dBxdx_ne*dx_ne+dBxdy_ne*dy_ne+(d2Bxdxdx_ne*dx2_ne+d2Bxdxdy_ne*dxdy_ne+d2Bxdydy_ne*dy2_ne);
    By_ne=By0_ne+dBydx_ne*dx_ne+dBydy_ne*dy_ne+(d2Bydxdx_ne*dx2_ne+d2Bydxdy_ne*dxdy_ne+d2Bydydy_ne*dy2_ne);
    
    *Bx= (dx_ne*dy_ne*Bx_sw-dx_nw*dy_nw*Bx_se-dx_se*dy_se*Bx_nw+dx_sw*dy_sw*Bx_ne)/(deltaxdeltay);
    *By= (dx_ne*dy_ne*By_sw-dx_nw*dy_nw*By_se-dx_se*dy_se*By_nw+dx_sw*dy_sw*By_ne)/(deltaxdeltay);
    
    /*printf("dBx_sw=%e, dBx_se=%e, dBx_nw=%e, dBx_ne=%e, Bxtot=%e \n",(Bx_sw-*Bx) / *Bx,(Bx_se-*Bx) / *Bx,(Bx_nw-*Bx) / *Bx,(Bx_ne-*Bx) / *Bx,*Bx);
    printf("dBy_sw=%e, dBy_se=%e, dBy_nw=%e, dBy_ne=%e, Bytot=%e \n",(By_sw-*By) / *By,(By_se-*By) / *By,(By_nw-*By) / *By,(By_ne-*By) / *By,*By);*/
    /*printf("(%e-%e-%e+%e)/%e=%f\n",dx_ne*dy_ne,dx_nw*dy_nw,dx_se*dy_se,dx_sw*dy_sw,deltaxdeltay, (dx_ne*dy_ne-dx_nw*dy_nw-dx_se*dy_se+dx_sw*dy_sw)/deltaxdeltay);*/
    /**Bx= Bx_sw;
    *By= By_sw;*/
    
}

void interpolate1ord(struct elem *Elem,int xi,int yi,double dx,double dy,double *Bx,double *By)
{
    int Nx=Elem->Nx;
    int Ny=Elem->Ny;
    double Bx0, dBxdx, dBxdy;
    double By0, dBydx, dBydy;
                
    Bx0=Elem->LUT_Bx[yi*Nx+xi];
    dBxdx=Elem->LUT_dBxdx[yi*Nx+xi];
    dBxdy=Elem->LUT_dBxdy[yi*Nx+xi];
    /*printf("dBxdx = %f \n", dBxdx);
    printf("dBxdy = %f \n", dBxdy);*/
    
    By0=Elem->LUT_By[yi*Nx+xi];
    dBydx=Elem->LUT_dBydx[yi*Nx+xi];
    dBydy=Elem->LUT_dBydy[yi*Nx+xi];
    
    *Bx=Bx0+dBxdx*dx+dBxdy*dy;
    *By=By0+dBydx*dx+dBydy*dy;
    /*printf("Bx = %f = %f+%f*%f+%f*%f \n", *Bx,Bx0,dBxdx,dx,dBxdy,dy);
    printf("By = %f = %f+%f*%f+%f*%f \n", *By,By0,dBydx,dx,dBydy,dy);*/
}

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


static void strthinkickFMrad(double* r, double ReSum, double ImSum, double L, double E0)
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
   /*double ReSum = B[max_order];
   double ImSum = A[max_order];
   double ReSumTemp;*/
   double irho=0;/*straight elements no curvature.*/
   double x ,xpr, y, ypr, p_norm,dp_0, B2P;
   double CRAD = CGAMMA*E0*E0*E0/(TWOPI*1e27);	/* [m]/[GeV^3] M.Sands (4.1) */
   
   /*for (i=max_order-1; i>=0; i--) {
      ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
      ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
      ReSum = ReSumTemp;
   }*/
   
   printf("ReSum= %f, ImSum= %f \n", ReSum, ImSum);
   
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
