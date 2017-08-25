#include <math.h>

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

/* from elegant code */
static void quadPartialFringeMatrix(double R[6][6], double K1, double inFringe, double *fringeInt, int part)
{
  double J1x, J2x, J3x, J1y, J2y, J3y;
  double K1sqr, expJ1x, expJ1y;

  R[4][4] = R[5][5] = 1;
  
  K1sqr = K1*K1;

  if (part==1) {
    J1x = inFringe*(K1*fringeInt[1] - 2*K1sqr*fringeInt[3]/3.);
    J2x = inFringe*(K1*fringeInt[2]);
    J3x = inFringe*(K1sqr*(fringeInt[2] + fringeInt[4]));

    K1 = -K1;
    J1y = inFringe*(K1*fringeInt[1] - 2*K1sqr*fringeInt[3]/3.);
    J2y = -J2x;
    J3y = J3x;
  } else {
    J1x = inFringe*(K1*fringeInt[1] + K1sqr*fringeInt[0]*fringeInt[2]/2);
    J2x = inFringe*(K1*fringeInt[2]);
    J3x = inFringe*(K1sqr*(fringeInt[4]-fringeInt[0]*fringeInt[1]));

    K1 = -K1;
    J1y = inFringe*(K1*fringeInt[1] + K1sqr*fringeInt[0]*fringeInt[2]);
    J2y = -J2x;
    J3y = J3x;
  }

  expJ1x = R[0][0] = exp(J1x);
  R[0][1] = J2x/expJ1x;
  R[1][0] = expJ1x*J3x;
  R[1][1] = (1 + J2x*J3x)/expJ1x;
  
  expJ1y = R[2][2] = exp(J1y);
  R[2][3] = J2y/expJ1y;
  R[3][2] = expJ1y*J3y;
  R[3][3] = (1 + J2y*J3y)/expJ1y;

  return;
}

static void linearQuadFringeElegantEntrance(double* r6, double b2, double *fringeIntM0, double *fringeIntP0)
{
    double R[6][6];
    double *fringeIntM, *fringeIntP;
    double delta, inFringe;
    /* quadrupole linear fringe field, from elegant code */
    inFringe=-1.0;
    fringeIntM = fringeIntP0;
    fringeIntP = fringeIntM0;
    delta = r6[4];
    /* determine first linear matrix for this delta */
    quadPartialFringeMatrix(R, b2/(1+delta), inFringe, fringeIntM, 1);
    r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
    r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
    r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
    r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
    /* nonlinear fringe field */
    QuadFringePassP(r6,b2);   /*This is original AT code*/
    /*Linear fringe fields from elegant*/
    inFringe=-1.0;
    /* determine and apply second linear matrix, from elegant code */
    quadPartialFringeMatrix(R, b2/(1+delta), inFringe, fringeIntP, 2);
    r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
    r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
    r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
    r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
}


static void linearQuadFringeElegantExit(double* r6, double b2, double *fringeIntM0, double *fringeIntP0)
{
    double R[6][6];
    double *fringeIntM, *fringeIntP;
    double delta, inFringe;
    /* quadrupole linear fringe field, from elegant code */
    inFringe=1.0;
    fringeIntM = fringeIntM0;
    fringeIntP = fringeIntP0;
    delta = r6[4];
    /* determine first linear matrix for this delta */
    quadPartialFringeMatrix(R, b2/(1+delta), inFringe, fringeIntM, 1);
    r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
    r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
    r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
    r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
    /* nonlinear fringe field */
    QuadFringePassN(r6,b2);   /*This is original AT code*/
    /*Linear fringe fields from elegant*/
    inFringe=1.0;
    /* determine and apply second linear matrix, from elegant code */
    quadPartialFringeMatrix(R, b2/(1+delta), inFringe, fringeIntP, 2);
    r6[0] = R[0][0]*r6[0] + R[0][1]*r6[1];
    r6[1] = R[1][0]*r6[0] + R[1][1]*r6[1];
    r6[2] = R[2][2]*r6[2] + R[2][3]*r6[3];
    r6[3] = R[3][2]*r6[2] + R[3][3]*r6[3];
}


