/*
 *----------------------------------------------------------------------------
 * Modification Log:
 * -----------------
 * .07  2024-05-09      J. Arenillas, ALBA, jarenillas@axt.email
 *                      New expression for the comptation of radiation kicks.
 *						Adding new functions for the computation of wiggler diffusion matrix.
 *						Extending GwigB to compute the z-coordinate of the magnetic field.
 *
 * .06  2018-06-01      A.Mash'al, ILSF, a-mashal@ipm.ir
 *                      Implementing Hessian of the Hamiltonian for computing wiggler transfer matrix
 *
 * .05  2018-02-13      O.Jimenez, ALBA Synchrotron, oscar.jimenez.1996@gmail.com
 *                      Implementing radiation loss.
 *
 * .04  2003-04-29      YK Wu, Duke University
 *			            Using scientific notation for constants.
 *                      Checked with TRACY pascal code.
 *                      Computing differential pathlength only.
 *
 * .03  2003-04-28      YK Wu, Duke University
 *			            Convert to C code and cross-checked with the pascal version;
 *
 * .02  2001-12-xx      Y. K. Wu, Duke University
 *                      Implementing DA version of the wiggler integrator for Pascal.
 *                      Gauge is disabled !!! (Dec. 4, 2001)
 *
 * .01  2001-02-12      Y. K. Wu, LBNL
 *                      Implementing a generic wiggler integrator
 *                      for paraxial-ray Hamiltonian approximation.
 *
 *
 *----------------------------------------------------------------------------
 *  Accelerator Physics Group, Duke FEL Lab, www.fel.duke.edu
 */

#ifndef  GWIG
#include "gwig.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#endif

/* struct used for GWigSymplecticRadPass */
struct gwigR{
  int Pmethod;      /* Integration Method */
  int PN;           /* Number of integration steps */
  double E0;        /* Energy of ring, [GeV] */
  double PB0;       /* B0 in [Tesla] */
  int Nw;           /* Number of periods */
  double Lw;        /* Wiggler Period [m] */
  int NHharm;       /* No. of horizontal harmonics */
  int NVharm;       /* No. of vertical harmonics */
  double Aw;        /* Wiggler parameter */
  double Zw;        /* Longitudinal variable [m] */
  double zStartH;
  double zStartV;  /* Start and end z coordinates of the wiggler field, which are computed */
  double zEndH;
  double zEndV;      /* based on the phase of the first harmonic to get matched dispersion. */
  double srCoef;
  double Po;        /* beta*gamma for reference particle */
  int HSplitPole;
  int VSplitPole;

  double HCw[WHmax];
  double VCw[WHmax];
  double HCw_raw[WHmax];
  double VCw_raw[WHmax];
  double Hkx[WHmax];
  double Hky[WHmax];
  double Hkz[WHmax];
  double Htz[WHmax];
  double Vkx[WHmax];
  double Vky[WHmax];
  double Vkz[WHmax];
  double Vtz[WHmax];
};

static const double q_e    = 1.602176462e-19; /* electron charge, [C] */
static const double m_e    = 9.10938188e-31;  /* electron mass, [kg] */
static const double clight = 2.99792458e8;    /* speed of light [m/s] */
static const double r_e    = 2.817940285e-15; /* electron classic radius,[m]*/
static const double XMC2   = 0.510998902e-03; /* mc^2 in GeV */
static const double PI     = 3.141592653589793238462643383279502884197e0;
static const double epsilon_o = 8.854187817e-12;  /*Vacuum permittivity*/

#define SQR(X) ((X)*(X))

static void GWigGauge(struct gwigR *pWig, double *X, int flag);
static void GWigPass_2nd(struct gwigR *pWig, double *X);
static void GWigPass_4th(struct gwigR *pWig, double *X);
static void GWigMap_2nd(struct gwigR *pWig, double *X, double dl);
static void GWigAx(struct gwigR *pWig, double *Xvec, double *pax, double *paxpy);
static void GWigAy(struct gwigR *pWig, double *Xvec, double *pay, double *paypx);
static void GWigRadiationKicks(struct gwigR *pWig, double *X, double *Bxyz, double dl);
static void GWigB(struct gwigR *pWig, double *Xvec, double *B);
static void AxHessian(struct gwigR *pWig, double *Xvec, double *pax);
static void AyHessian(struct gwigR *pWig, double *Xvec, double *pay);
static void Hessian(struct gwigR *pWig, double *Xvec, double *H2);
static void GWigInit2(struct gwigR *pWig, double design_energy, double Ltot, double Lw, double Bmax,
	      int Nstep, int Nmeth, int NHharm, int NVharm,
		  int HSplitPole, int VSplitPole, double *zEndPointH, double *zEndPointV,
	      double *By, double *Bx, double *T1, double *T2,
	      double *R1, double *R2);
static double sinc(double x);


static void GWigGauge(struct gwigR *pWig, double *X, int flag)
{
  double ax, ay, axpy, aypx;
  GWigAx(pWig, X, &ax, &axpy);
  GWigAy(pWig, X, &ay, &aypx);

  if (flag == Elem_Entrance) {
    /* At the entrance of the wiggler */
    X[1] = X[1] + ax;
    X[3] = X[3] + ay;
  } else if (flag == Elem_Exit) {
    /* At the exit of the wiggler */
    X[1] = X[1] - ax;
    X[3] = X[3] - ay;
  } else {
    printf("  GWigGauge: Unknown flag = %i\n", flag);
  }
}


static void GWigPass_2nd(struct gwigR *pWig, double *X)
{
  int    i, Nstep;
  double dl;
  double B[3];
  double ax, ay, axpy, aypx;

  Nstep = pWig->PN*(pWig->Nw);
  dl    = pWig->Lw/(pWig->PN);

  GWigAx(pWig, X, &ax, &axpy);
  GWigAy(pWig, X, &ay, &aypx);
  GWigB(pWig, X, B);
  X[1] -= ax;
  X[3] -= ay;
  GWigRadiationKicks(pWig, X, B, dl);
  X[1] += ax;
  X[3] += ay;

  for (i = 1; i <= Nstep; i++) {
    GWigMap_2nd(pWig, X, dl);
    GWigAx(pWig, X, &ax, &axpy);
    GWigAy(pWig, X, &ay, &aypx);
    GWigB(pWig, X, B);
    X[1] -= ax;
    X[3] -= ay;
	  GWigRadiationKicks(pWig, X, B, dl);
	  X[1] += ax;
    X[3] += ay;
  }
}


static void GWigPass_4th(struct gwigR *pWig, double *X)
{
  const double x1 = 1.3512071919596576340476878089715e0;
  const double x0 =-1.7024143839193152680953756179429e0;

  int    i, Nstep;
  double dl, dl1, dl0;
  double B[3];
  double ax, ay, axpy, aypx;

  Nstep = pWig->PN*(pWig->Nw);
  dl = pWig->Lw/(pWig->PN);

  dl1 = x1*dl;
  dl0 = x0*dl;

  GWigAx(pWig, X, &ax, &axpy);
  GWigAy(pWig, X, &ay, &aypx);
  GWigB(pWig, X, B);
  X[1] -= ax;
  X[3] -= ay;
  GWigRadiationKicks(pWig, X, B, dl);
  X[1] += ax;
  X[3] += ay;

  for (i = 1; i <= Nstep; i++) {
    GWigMap_2nd(pWig, X, dl1);
    GWigMap_2nd(pWig, X, dl0);
    GWigMap_2nd(pWig, X, dl1);
    GWigAx(pWig, X, &ax, &axpy);
    GWigAy(pWig, X, &ay, &aypx);
    GWigB(pWig, X, B);
    X[1] -= ax;
    X[3] -= ay;
    GWigRadiationKicks(pWig, X, B, dl);
    X[1] += ax;
    X[3] += ay;
  }
}


static void GWigMap_2nd(struct gwigR *pWig, double *X, double dl)
{
  double dld, dl2, dl2d;
  double ax, ay, axpy, aypx;

  dld  = dl/(1.0e0 + X[4]);
  dl2  = 0.5e0 * dl;
  dl2d = dl2/(1.0e0 + X[4]);

  /* Step1: increase a half step in z */
  pWig->Zw = pWig->Zw + dl2;

  /* Step2: a half drift in y */
  GWigAy(pWig, X, &ay, &aypx);
  X[1] = X[1] - aypx;
  X[3] = X[3] - ay;

  X[2] = X[2] + dl2d*X[3];
  X[5] = X[5] + 0.5e0*dl2d*(X[3]*X[3])/(1.0e0+X[4]);

  GWigAy(pWig, X, &ay, &aypx);
  X[1] = X[1] + aypx;
  X[3] = X[3] + ay;

  /* Step3: a full drift in x */
  GWigAx(pWig, X, &ax, &axpy);
  X[1] = X[1] - ax;
  X[3] = X[3] - axpy;

  X[0] = X[0] + dld*X[1];
/* Full path length
  X[5] = X[5] + dl + 0.5e0*dld*(X[1]*X[1])/(1.0e0+X[4]);
*/
  /* Differential path length only */
  X[5] = X[5] + 0.5e0*dld*(X[1]*X[1])/(1.0e0+X[4]);

  GWigAx(pWig, X, &ax, &axpy);
  X[1] = X[1] + ax;
  X[3] = X[3] + axpy;

  /* Step4: a half drift in y */
  GWigAy(pWig, X, &ay, &aypx);
  X[1] = X[1] - aypx;
  X[3] = X[3] - ay;

  X[2] = X[2] + dl2d*X[3];
  X[5] = X[5] + 0.5e0*dl2d*(X[3]*X[3])/(1.0e0+X[4]);

  GWigAy(pWig, X, &ay, &aypx);
  X[1] = X[1] + aypx;
  X[3] = X[3] + ay;

  /* Step5: increase a half step in z */
  pWig->Zw = pWig->Zw + dl2;
}


static void GWigAx(struct gwigR *pWig, double *Xvec, double *pax, double *paxpy)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sxkx, chx, shx;
  double cy, sy, chy, shy, sz;
  double gamma0, beta0;
  double ax, axpy;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = 2e0*PI/(pWig->Lw);
  ax   = 0e0;
  axpy = 0e0;

  gamma0   = pWig->E0/XMC2;
  beta0    = sqrt(1e0 - 1e0/(gamma0*gamma0));
  pWig->Aw = (q_e/m_e/clight)/(2e0*PI) * (pWig->Lw) * (pWig->PB0);

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for (i = 0; i < pWig->NHharm; i++) {
    pWig->HCw[i] = pWig->HCw_raw[i]*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    cx  = cos(kx*x);
    chy = cosh(ky * y);
    sz  = sin(kz * z + tz);
    ax  = ax + (pWig->HCw[i])*(kw/kz)*cx*chy*sz;

    shy = sinh(ky * y);
    if (fabs(kx/kw) > GWIG_EPS) {
      sxkx = sin(kx * x)/kx;
    } else {
      sxkx = x*sinc(kx*x);
    }

    axpy = axpy + pWig->HCw[i]*(kw/kz)*ky*sxkx*shy*sz;
  }


  /* Vertical Wiggler: note that one potentially could have: ky=0 */
  for (i = 0; i < pWig->NVharm; i++) {
    pWig->VCw[i] = pWig->VCw_raw[i]*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Vkx[i];
    ky = pWig->Vky[i];
    kz = pWig->Vkz[i];
    tz = pWig->Vtz[i];

    shx = sinh(kx * x);
    sy  = sin(ky * y);
    sz  = sin(kz * z + tz);
    ax  = ax + pWig->VCw[i]*(kw/kz)*(ky/kx)*shx*sy*sz;

    chx = cosh(kx * x);
    cy  = cos(ky * y);
    axpy = axpy + pWig->VCw[i]*(kw/kz)* pow(ky/kx,2) *chx*cy*sz;
  }

  *pax   = ax;
  *paxpy = axpy;
}


static void GWigAy(struct gwigR *pWig, double *Xvec, double *pay, double *paypx)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sx, chx, shx;
  double cy, syky, chy, shy, sz;
  double gamma0, beta0;
  double ay, aypx;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = 2e0*PI/(pWig->Lw);
  ay   = 0e0;
  aypx = 0e0;

  gamma0  = pWig->E0/XMC2;
  beta0   = sqrt(1e0 - 1e0/(gamma0*gamma0));
  pWig->Aw = (q_e/m_e/clight)/(2e0*PI) * (pWig->Lw) * (pWig->PB0);

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for (i = 0; i < pWig->NHharm; i++){
    pWig->HCw[i] = (pWig->HCw_raw[i])*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    sx = sin(kx * x);
    shy = sinh(ky * y);
    sz  = sin(kz * z + tz);
    ay  = ay + (pWig->HCw[i])*(kw/kz)*(kx/ky)*sx*shy*sz;

    cx  = cos(kx * x);
    chy = cosh(ky * y);

    aypx = aypx + (pWig->HCw[i])*(kw/kz)*pow(kx/ky,2) * cx*chy*sz;
  }

  /* Vertical Wiggler: note that one potentially could have: ky=0 */
  for (i = 0; i < pWig->NVharm; i++) {
    pWig->VCw[i] = (pWig->VCw_raw[i])*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Vkx[i];
    ky = pWig->Vky[i];
    kz = pWig->Vkz[i];
    tz = pWig->Vtz[i];

    chx = cosh(kx * x);
    cy  = cos(ky * y);
    sz  = sin(kz * z + tz);
    ay  = ay + (pWig->VCw[i])*(kw/kz)*chx*cy*sz;

    shx = sinh(kx * x);
    if (fabs(ky/kw) > GWIG_EPS) {
      syky  = sin(ky * y)/ky;
    } else {
      syky = y * sinc(ky * y);
    }
    aypx = aypx + (pWig->VCw[i])*(kw/kz)* kx*shx*syky*sz;
   }

  *pay = ay;
  *paypx = aypx;
}


static double sinc(double x)
{
  double x2, result;
/* Expand sinc(x) = sin(x)/x to x^8 */
  x2 = x*x;
  result = 1e0 - x2/6e0*(1e0 - x2/20e0 *(1e0 - x2/42e0*(1e0-x2/72e0)));
  return result;
}


static void GWigInit2(struct gwigR *Wig,double design_energy, double Ltot, double Lw,
            double Bmax, int Nstep, int Nmeth, int NHharm, int NVharm,
            int HSplitPole, int VSplitPole, double *zEndPointH,
            double *zEndPointV, double *By, double *Bx, double *T1,
            double *T2, double *R1, double *R2)
{
    double *tmppr;
    int    i;
    double kw;

    Wig->E0 = design_energy;
    Wig->Po = Wig->E0/XMC2;
    Wig->Pmethod = Nmeth;
    Wig->PN = Nstep;
    Wig->Nw = (int)(Ltot / Lw);
    Wig->NHharm = NHharm;
    Wig->NVharm = NVharm;
    Wig->PB0 = Bmax;
    Wig->Lw = Lw;

    /*------------------ radiation including -------------------*/
    Wig->srCoef = (q_e*q_e)*((Wig->Po)*(Wig->Po)*(Wig->Po))/(6*PI*epsilon_o*m_e*(clight*clight));
    Wig->HSplitPole = HSplitPole;
    Wig->VSplitPole = VSplitPole;
    Wig->zStartH = zEndPointH[0];
    Wig->zEndH = zEndPointH[1];
    Wig->zStartV = zEndPointV[0];
    Wig->zEndV = zEndPointV[1];
    /*----------------------------------------------------------*/

    kw = 2.0e0*PI/(Wig->Lw);
    Wig->Zw = 0.0;
    Wig->Aw = 0.0;
    tmppr = By;
    for (i = 0; i < NHharm; i++) {
        tmppr++;
        Wig->HCw[i] = 0.0;
        Wig->HCw_raw[i] = *tmppr;
        tmppr++;
        Wig->Hkx[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Hky[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Hkz[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Htz[i] = *tmppr;
        tmppr++;
    }
    tmppr = Bx;
    for (i = 0; i < NVharm; i++) {
        tmppr++;
        Wig->VCw[i] = 0.0;
        Wig->VCw_raw[i] = *tmppr;
        tmppr++;
        Wig->Vkx[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Vky[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Vkz[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Vtz[i] = *tmppr;
        tmppr++;
    }
    for (i = NHharm; i< WHmax; i++) {
        Wig->HCw[i] = 0.0;
        Wig->HCw_raw[i] = 0.0;
        Wig->Hkx[i] = 0.0;
        Wig->Hky[i] = 0.0;
        Wig->Hkz[i] = 0.0;
        Wig->Htz[i] = 0.0;
    }
    for (i = NVharm; i< WHmax; i++) {
        Wig->VCw[i] = 0.0;
        Wig->VCw_raw[i] = 0.0;
        Wig->Vkx[i] = 0.0;
        Wig->Vky[i] = 0.0;
        Wig->Vkz[i] = 0.0;
        Wig->Vtz[i] = 0.0;
    }
}


static void GWigB(struct gwigR *pWig, double *Xvec, double *B)
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


static void GWigRadiationKicks(struct gwigR *pWig, double *X, double *Bxy, double dl)
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
  H = (pWig->Po)/586.679074042074490;

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


/* Hessian functions */
static void Hessian(struct gwigR *pWig, double *Xvec, double *H2)
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


static void AxHessian(struct gwigR *pWig, double *Xvec, double *pax)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, chx, shx;
  double sx;
  double cy, sy, chy, shy, sz;
  double gamma0, beta0;
  double ax,axx,axxx,axy,axyy,axxy;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = 2e0*PI/(pWig->Lw);
  ax   = 0e0;
  axx  = 0e0;
  axxx = 0e0;
  axy  = 0e0;
  axyy = 0e0;
  axxy = 0e0;
  gamma0   = pWig->Po;
  beta0    = sqrt(1e0 - 1e0/(gamma0*gamma0));
  pWig->Aw = (q_e/m_e/clight)/(2e0*PI) * (pWig->Lw) * (pWig->PB0);

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for (i = 0; i < pWig->NHharm; i++) {
    pWig->HCw[i] = pWig->HCw_raw[i]*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    cx  = cos(kx*x);
	sx  = sin(kx*x);
    chy = cosh(ky * y);
    shy = sinh(ky * y);
    sz  = sin(kz * z + tz);
    ax  = ax + (pWig->HCw[i])*(kw/kz)*cx*chy*sz;
	axx = axx- (pWig->HCw[i])*(kx)*(kw/kz)*sx*chy*sz;
	axxx=axxx- (pWig->HCw[i])*(kx*kx)*(kw/kz)*cx*chy*sz;
	axy = axy+ (pWig->HCw[i])*(ky)*(kw/kz)*cx*shy*sz;
	axyy=axyy+ (pWig->HCw[i])*(ky*ky)*(kw/kz)*cx*chy*sz;
	axxy=axxy- (pWig->HCw[i])*(kx*ky)*(kw/kz)*sx*shy*sz;
  }

  /* Vertical Wiggler: note that one potentially could have: ky=0 */
  for (i = 0; i < pWig->NVharm; i++) {
    pWig->VCw[i] = pWig->VCw_raw[i]*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Vkx[i];
    ky = pWig->Vky[i];
    kz = pWig->Vkz[i];
    tz = pWig->Vtz[i];

    shx = sinh(kx * x);
    sy  = sin(ky * y);
	chx = cosh(kx * x);
    cy  = cos(ky * y);
    sz  = sin(kz * z + tz);
    ax  = ax + pWig->VCw[i]*(kw/kz)*(ky/kx)*shx*sy*sz;
	axx = axx+ pWig->VCw[i]*(kw/kz)*(ky)*chx*sy*sz;
	axxx=axxx+ pWig->VCw[i]*(kw/kz)*(ky*kx)*chx*sy*sz;
	axy = axy+ pWig->VCw[i]*(kw/kz)*(ky)*(ky/kx)*shx*cy*sz;
	axyy=axyy- pWig->VCw[i]*(kw/kz)*(ky*ky)*(ky/kx)*shx*sy*sz;
	axxy=axxy+ pWig->VCw[i]*(kw/kz)*(ky*ky)*chx*cy*sz;
  }

  pax[0]   = ax;
  pax[1]   = axx;
  pax[2]   = axxx;
  pax[3]   = axy;
  pax[4]   = axyy;
  pax[5]   = axxy;
}


static void AyHessian(struct gwigR *pWig, double *Xvec, double *pay)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sx, chx, shx, sy;
  double cy, chy, shy, sz;
  double gamma0, beta0;
  double ay, ayx, ayxx, ayy, ayyy, ayxy;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = 2e0*PI/(pWig->Lw);
  ay   = 0e0;
  ayx  = 0e0;
  ayxx = 0e0;
  ayy  = 0e0;
  ayyy = 0e0;
  ayxy = 0e0;

  gamma0  = pWig->Po;
  beta0   = sqrt(1e0 - 1e0/(gamma0*gamma0));
  pWig->Aw = (q_e/m_e/clight)/(2e0*PI) * (pWig->Lw) * (pWig->PB0);

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for (i = 0; i < pWig->NHharm; i++){
    pWig->HCw[i] = (pWig->HCw_raw[i])*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    sx = sin(kx * x);
    shy = sinh(ky * y);
    sz  = sin(kz * z + tz);
	cx  = cos(kx * x);
    chy = cosh(ky * y);


    ay  = ay + (pWig->HCw[i])*(kw/kz)*(kx/ky)*sx*shy*sz;
	ayx = ayx+ (pWig->HCw[i])*(kx)*(kw/kz)*(kx/ky)*cx*shy*sz;
    ayxx=ayxx- (pWig->HCw[i])*(kw/kz)*(kx*kx)*(kx/ky)*sx*shy*sz;
	ayy = ayy+ (pWig->HCw[i])*(kw/kz)*(kx)*sx*chy*sz;
	ayyy=ayyy+ (pWig->HCw[i])*(kw/kz)*(ky*kx)*sx*shy*sz;
	ayxy=ayxy+ (pWig->HCw[i])*(kw/kz)*(kx*kx)*cx*chy*sz;
  }

  /* Vertical Wiggler: note that one potentially could have: ky=0 */
  for (i = 0; i < pWig->NVharm; i++) {
    pWig->VCw[i] = (pWig->VCw_raw[i])*(pWig->Aw)/(gamma0*beta0);
    kx = pWig->Vkx[i];
    ky = pWig->Vky[i];
    kz = pWig->Vkz[i];
    tz = pWig->Vtz[i];

    chx = cosh(kx * x);
    cy  = cos(ky * y);
	sy  = sin(ky * y);
    sz  = sin(kz * z + tz);
	shx = sinh(kx * x);

    ay  = ay + (pWig->VCw[i])*(kw/kz)*chx*cy*sz;
	ayx = ayx+ (pWig->VCw[i])*(kx)*(kw/kz)*shx*cy*sz;
	ayxx=ayxx+ (pWig->VCw[i])*(kx*kx)*(kw/kz)*chx*cy*sz;
	ayy = ayy- (pWig->VCw[i])*(ky)*(kw/kz)*chx*sy*sz;
	ayyy=ayyy- (pWig->VCw[i])*(ky*ky)*(kw/kz)*chx*cy*sz;
	ayxy=ayxy- (pWig->VCw[i])*(kx*ky)*(kw/kz)*shx*sy*sz;
  }

  pay[0]   = ay;
  pay[1]   = ayx;
  pay[2]   = ayxx;
  pay[3]   = ayy;
  pay[4]   = ayyy;
  pay[5]   = ayxy;
}

static void print66(const char *title, double *m)
{
    atPrintf("\n%s\n", title);
    for (int i=0; i<6; i++) {
        int k=6*i;
        atPrintf("%.12g, %.12g, %.12g, %.12g, %.12g, %.12g\n", m[k+0], m[k+1], m[k+2], m[k+3], m[k+4], m[k+5]);
    }
}
