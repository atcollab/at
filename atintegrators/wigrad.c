#include "gwig.c"

static void GWigB(struct gwig *pWig, double *Xvec, double *B)
/* Compute magnetic field at particle location.
 * Added by M. Borland, August 2007.
 */
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz;
  double cx, sx, chx, shx;
  double cy, sy, chy, shy;
  double cz, sz;
  /* B0 is a reserved symbol on MacOS, defined in termios.h */
  double _B0;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  B[0] = 0;
  B[1] = 0;
  B[2] = 0;

  if (pWig->NHharm && z>=pWig->zStartH && z<=pWig->zEndH) {
    _B0 = pWig->PB0;
    if (!pWig->HSplitPole) {
      /* Normal Horizontal Wiggler: note that one potentially could have: kx=0 */
      for (i = 0; i < pWig->NHharm; i++) {
        kx = pWig->Hkx[i];
        ky = pWig->Hky[i];
        kz = pWig->Hkz[i];
        tz = pWig->Htz[i];

        sx  = sin(kx*x);
        cx  = cos(kx*x);
        chy = cosh(ky * y);
        shy = sinh(ky * y);
        cz = cos(kz*z+tz);
		sz = sin(kz*z+tz);

        /* Accumulate field values in user-supplied array (Bx, By, Bz) */
        B[0] += _B0*pWig->HCw_raw[i]*kx/ky*sx*shy*cz;
        B[1] -= _B0*pWig->HCw_raw[i]*cx*chy*cz;
		B[2] += _B0*pWig->HCw_raw[i]*kz/ky*cx*shy*sz;
      }
    } else {
      /* Split-pole Horizontal Wiggler: note that one potentially could have: ky=0 (caught in main routine) */
      for (i = 0; i < pWig->NHharm; i++) {
        kx = pWig->Hkx[i];
        ky = pWig->Hky[i];
        kz = pWig->Hkz[i];
        tz = pWig->Htz[i];

        shx = sinh(kx*x);
        chx = cosh(kx*x);
        cy  = cos(ky * y);
        sy  = sin(ky * y);
        cz  = cos(kz*z+tz);
		cx  = cos(kx*x);
		chy = cosh(ky * y);
        sz = sin(kz*z+tz);

		/* Accumulate field values in user-supplied array (Bx, By, Bz) */
        B[0] -= _B0*pWig->HCw_raw[i]*kx/ky*shx*sy*cz;
        B[1] -= _B0*pWig->HCw_raw[i]*chx*cy*cz;
		B[2] -= _B0*pWig->HCw_raw[i]*kz/ky*cx*chy*sz;
      }
    }
  }

  if (pWig->NVharm && z>=pWig->zStartV && z<=pWig->zEndV) {
    _B0 = pWig->PB0;
    if (!pWig->VSplitPole) {
      /* Normal Vertical Wiggler: note that one potentially could have: ky=0 */
      for (i = 0; i < pWig->NVharm; i++) {
        kx = pWig->Vkx[i];
        ky = pWig->Vky[i];
        kz = pWig->Vkz[i];
        tz = pWig->Vtz[i];

        shx = sinh(kx * x);
        chx = cosh(kx * x);
        sy  = sin(ky * y);
        cy  = cos(ky * y);
        cz  = cos(kz*z + tz);
		sz = sin(kz*z+tz);

        /* Accumulate field values in user-supplied array (Bx, By, Bz) */
        B[0] += _B0*pWig->VCw_raw[i]*chx*cy*cz;
        B[1] -= _B0*pWig->VCw_raw[i]*ky/kx*shx*sy*cz;
		B[2] -= _B0*pWig->VCw_raw[i]*kz/kx*cy*shx*sz;
      }
    } else {
      /* Split-pole Vertical Wiggler: note that one potentially could have: kx=0 (caught in main routine) */
      for (i = 0; i < pWig->NVharm; i++) {
        kx = pWig->Vkx[i];
        ky = pWig->Vky[i];
        kz = pWig->Vkz[i];
        tz = pWig->Vtz[i];

        sx  = sin(kx * x);
        cx  = cos(kx * x);
        shy = sinh(ky * y);
        chy = cosh(ky * y);
        cz  = cos(kz*z + tz);

        /* Accumulate field values in user-supplied array (Bx, By, Bz) */
        B[0] += _B0*pWig->VCw_raw[i]*cx*chy*cz;
        B[1] += _B0*pWig->VCw_raw[i]*ky/kx*sx*shy*cz;
		B[2] += _B0*pWig->VCw_raw[i]*kz/kx*cx*shy*sz;
      }
    }
  }
}

static void GWigRadiationKicks(struct gwig *pWig, double *X, double *Bxy, double dl)
/* Apply kicks for synchrotron radiation.
 * Added by M. Borland, August 2007.
 */
{
  double irho2, H, dFactor;
  double B2;
  double dDelta;

  /* B^2 in T^2 */
  B2 = (Bxy[0]*Bxy[0]) + (Bxy[1]*Bxy[1]);
  if (B2==0)
	  return;

  /* Beam rigidity in T*m */
  H = 1.0e9 * pWig->Po * __E0 / C0;

  /* 1/rho^2 */
  irho2 = B2/(H*H);

  /* (1+delta)^2 */
  dFactor = ((1+X[4])*(1+X[4]));

  /* Classical radiation loss */
    dDelta = -(pWig->srCoef)*dFactor*irho2*dl;
    X[4] += dDelta;
    X[1] *= (1+dDelta);
    X[3] *= (1+dDelta);

}

static void AxHessian(struct gwig *pWig, double *Xvec, double *pax)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, chx, shx;
  double sx;
  double cy, sy, chy, shy, sz;
  double ax,axx,axxx,axy,axyy,axxy;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = TWOPI/(pWig->Lw);
  ax   = 0e0;
  axx  = 0e0;
  axxx = 0e0;
  axy  = 0e0;
  axyy = 0e0;
  axxy = 0e0;

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for (i = 0; i < pWig->NHharm; i++) {
    double HCw = pWig->HCw_raw[i] * pWig->Aw / pWig->Po;
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    cx  = cos(kx*x);
	sx  = sin(kx*x);
    chy = cosh(ky * y);
    shy = sinh(ky * y);
    sz  = sin(kz * z + tz);
    ax  = ax + HCw*(kw/kz)*cx*chy*sz;
	axx = axx- HCw*(kx)*(kw/kz)*sx*chy*sz;
	axxx=axxx- HCw*(kx*kx)*(kw/kz)*cx*chy*sz;
	axy = axy+ HCw*(ky)*(kw/kz)*cx*shy*sz;
	axyy=axyy+ HCw*(ky*ky)*(kw/kz)*cx*chy*sz;
	axxy=axxy- HCw*(kx*ky)*(kw/kz)*sx*shy*sz;
  }

  /* Vertical Wiggler: note that one potentially could have: ky=0 */
  for (i = 0; i < pWig->NVharm; i++) {
    double VCw = pWig->VCw_raw[i] * pWig->Aw / pWig->Po;
    kx = pWig->Vkx[i];
    ky = pWig->Vky[i];
    kz = pWig->Vkz[i];
    tz = pWig->Vtz[i];

    shx = sinh(kx * x);
    sy  = sin(ky * y);
	chx = cosh(kx * x);
    cy  = cos(ky * y);
    sz  = sin(kz * z + tz);
    ax  = ax + VCw*(kw/kz)*(ky/kx)*shx*sy*sz;
	axx = axx+ VCw*(kw/kz)*(ky)*chx*sy*sz;
	axxx=axxx+ VCw*(kw/kz)*(ky*kx)*chx*sy*sz;
	axy = axy+ VCw*(kw/kz)*(ky)*(ky/kx)*shx*cy*sz;
	axyy=axyy- VCw*(kw/kz)*(ky*ky)*(ky/kx)*shx*sy*sz;
	axxy=axxy+ VCw*(kw/kz)*(ky*ky)*chx*cy*sz;
  }

  pax[0]   = ax;
  pax[1]   = axx;
  pax[2]   = axxx;
  pax[3]   = axy;
  pax[4]   = axyy;
  pax[5]   = axxy;
}

static void AyHessian(struct gwig *pWig, double *Xvec, double *pay)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sx, chx, shx, sy;
  double cy, chy, shy, sz;
  double ay, ayx, ayxx, ayy, ayyy, ayxy;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = TWOPI/(pWig->Lw);
  ay   = 0e0;
  ayx  = 0e0;
  ayxx = 0e0;
  ayy  = 0e0;
  ayyy = 0e0;
  ayxy = 0e0;

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for (i = 0; i < pWig->NHharm; i++) {
    double HCw = pWig->HCw_raw[i] * pWig->Aw / pWig->Po;
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    sx = sin(kx * x);
    shy = sinh(ky * y);
    sz  = sin(kz * z + tz);
	cx  = cos(kx * x);
    chy = cosh(ky * y);


    ay  = ay + HCw*(kw/kz)*(kx/ky)*sx*shy*sz;
	ayx = ayx+ HCw*(kx)*(kw/kz)*(kx/ky)*cx*shy*sz;
    ayxx=ayxx- HCw*(kw/kz)*(kx*kx)*(kx/ky)*sx*shy*sz;
	ayy = ayy+ HCw*(kw/kz)*(kx)*sx*chy*sz;
	ayyy=ayyy+ HCw*(kw/kz)*(ky*kx)*sx*shy*sz;
	ayxy=ayxy+ HCw*(kw/kz)*(kx*kx)*cx*chy*sz;
  }

  /* Vertical Wiggler: note that one potentially could have: ky=0 */
  for (i = 0; i < pWig->NVharm; i++) {
    double VCw = pWig->VCw_raw[i] * pWig->Aw / pWig->Po;
    kx = pWig->Vkx[i];
    ky = pWig->Vky[i];
    kz = pWig->Vkz[i];
    tz = pWig->Vtz[i];

    chx = cosh(kx * x);
    cy  = cos(ky * y);
	sy  = sin(ky * y);
    sz  = sin(kz * z + tz);
	shx = sinh(kx * x);

    ay  = ay + VCw*(kw/kz)*chx*cy*sz;
	ayx = ayx+ VCw*(kx)*(kw/kz)*shx*cy*sz;
	ayxx=ayxx+ VCw*(kx*kx)*(kw/kz)*chx*cy*sz;
	ayy = ayy- VCw*(ky)*(kw/kz)*chx*sy*sz;
	ayyy=ayyy- VCw*(ky*ky)*(kw/kz)*chx*cy*sz;
	ayxy=ayxy- VCw*(kx*ky)*(kw/kz)*shx*sy*sz;
  }

  pay[0]   = ay;
  pay[1]   = ayx;
  pay[2]   = ayxx;
  pay[3]   = ayy;
  pay[4]   = ayyy;
  pay[5]   = ayxy;
}

/* Hessian functions */
static void Hessian(struct gwig *pWig, double *Xvec, double *H2)
{
	int i,j,k;
	double px, py,D;
	double ax,axx,axxx,axy,axxy;
	double ay, ayx, ayxx, ayy, ayxy;
	double Pax[6];
	double Pay[6];
	double H[6][6];

	AxHessian(pWig, Xvec, Pax);
	AyHessian(pWig, Xvec, Pay);

	ax  =Pax[0];
	axx =Pax[1];
	axxx=Pax[2];
	axy =Pax[3];
	axxy=Pax[5];

	ay  =Pay[0];
	ayx =Pay[1];
	ayxx=Pay[2];
	ayy =Pay[3];
	ayxy=Pay[5];

	px=Xvec[1];
	py=Xvec[3];
	D =Xvec[4];

	for(i=0;i<6;i++){
		for(j=0;j<6;j++){
			H[i][j]=0;
		}
	}

	H[0][0]= 0.5*(1/(1+D))*(((ax-px)*axxx)+ ((ay-py)*ayxx)+ (axx*axx)+ (ayx*ayx));
	H[0][1]=-0.5*(axx/(1+D));
	H[0][2]= 0.5*(1/(1+D))*(((ax-px)*axxy)+ ((ay-py)*ayxy)+ (axx*axy)+ (ayx*ayy));
	H[0][3]= -0.5*(ayx/(1+D));
	H[0][4]= -0.5*(1/((1+D)*(1+D)))*( ((ax-px)*axx)+ ((ay-py)*ayx));
	H[0][5]=  0.00;

	H[1][0]=-0.5*(axx/(1+D));
	H[1][1]= 0;
	H[1][2]=-0.5*(axy/(1+D));
	H[1][3]= 0;
	H[1][4]= 0.5*ax/((1+D)*(1+D));
	H[1][5]= 0.00;

	H[2][0]= 0.5*(1/(1+D))*(((ax-px)*axxy)+ ((ay-py)*ayxy)+ (axx*axy)+ (ayx*ayy));
	H[2][1]=-0.5*(axy/(1+D));
	H[2][2]= 0.5*(1/(1+D))*(((ax-px)*axxy)+ ((ay-py)*ayxy)+ (axy*axy)+ (ayy*ayy));
	H[2][3]= -0.5*(ayy/(1+D));
	H[2][4]= -0.5*(1/((1+D)*(1+D)))*( ((ax-px)*axy)+ ((ay-py)*ayy));
	H[2][5]=0;

	H[3][0]=-0.5*(ayx/(1+D));
	H[3][1]= 0;
	H[3][2]=-0.5*(ayy/(1+D));
	H[3][3]= 0;
	H[3][4]=0.5 *ay/((1+D)*(1+D));
	H[3][5]=0;

	H[4][0]= -0.5*(1/((1+D)*(1+D)))*( ((ax-px)*axx)+ ((ay-py)*ayx));
	H[4][1]= 0.5 *ax/((1+D)*(1+D));
	H[4][2]= -0.5*(1/((1+D)*(1+D)))*( ((ax-px)*axy)+ ((ay-py)*ayy));
	H[4][3]= 0.5 *ay/((1+D)*(1+D));
	H[4][4]=0.5 *(1/((1+D)*(1+D) *(1+D)))*(((ax*ax)+(ay*ay))-2*((ax*px)+(ay*py)));
	H[4][5]=0;

	for (i=0;i<6;i++){
		for(j=0;j<6;j++){
			k=j+i*6;
			H2[k]=H[i][j];
		}
	}
}
