
static void QuadFringePassP(double* r, const double b2)
{
/*	x=r[0],px=r[1],y=r[2],py=r[3],delta=r[4],ct=r[5] 
	Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam Dynamics..."by E. Forest */
   register double u = b2/(12.0*(1.0+r[4]));
   register double x2 = r[0]*r[0];
   register double z2 = r[2]*r[2];
   register double xz = r[0]*r[2];
   register double gx = u * (x2+3*z2) * r[0];
   register double gz = u * (z2+3*x2) * r[2];
   register double r1tmp=0;
   register double r3tmp=0;

   r[0]+=gx;
   r1tmp=3*u*(2*xz*r[3]-(x2+z2)*r[1]);
   
   r[2]-=gz;
   
   r3tmp=3*u*(2*xz*r[1]-(x2+z2)*r[3]);
   r[5]-=(gz*r[3] - gx*r[1])/(1+r[4]);
   
   r[1]+=r1tmp;
   r[3]-=r3tmp;
}

static void QuadFringePassN(double* r, const double b2)
{
/*	x=r[0],px=r[1],y=r[2],py=r[3],delta=r[4],ct=r[5] 
	Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam Dynamics..."by E. Forest */
   register double u = b2/(12.0*(1.0+r[4]));
   register double x2 = r[0]*r[0];
   register double z2 = r[2]*r[2];
   register double xz = r[0]*r[2];
   register double gx = u * (x2+3*z2) * r[0];
   register double gz = u * (z2+3*x2) * r[2];
   register double r1tmp=0;
   register double r3tmp=0;

   r[0]-=gx;
   r1tmp=3*u*(2*xz*r[3]-(x2+z2)*r[1]);
   
   r[2]+=gz;
   
   r3tmp=3*u*(2*xz*r[1]-(x2+z2)*r[3]);
   r[5]+=(gz*r[3] - gx*r[1])/(1+r[4]);
   
   r[1]-=r1tmp;
   r[3]+=r3tmp;
}

