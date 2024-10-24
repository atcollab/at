/*
 *----------------------------------------------------------------------------
 * Modification Log:
 * -----------------
 * .04  2003-04-29      YK Wu, Duke University
 *			using scientific notation for constants. 
 *                      Checked with TRACY pascal code.
 *                      Computing differential pathlength only.
 *
 * .03  2003-04-28      YK Wu, Duke University 
 *			Convert to C code and cross-checked with the pascal version;
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

#include "gwig.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

static void GWigInit2(struct gwig *Wig,double gamma, double Ltot, double Lw,
            double Bmax, int Nstep, int Nmeth, int NHharm, int NVharm,
            int HSplitPole, int VSplitPole, double *By, double *Bx,
            double *T1, double *T2, double *R1, double *R2)
{
    double *tmppr;
    int    i;
    double kw;

    Wig->Po = sqrt(gamma*gamma - 1.0);
    Wig->Pmethod = Nmeth;
    Wig->PN = Nstep;
    Wig->Nw = (int)(Ltot / Lw);
    Wig->NHharm = NHharm;
    Wig->NVharm = NVharm;
    Wig->PB0 = Bmax;
    Wig->Lw = Lw;

    /*------------------ radiation including -------------------*/
    Wig->srCoef = 2.0/3.0*__RE*pow(gamma, 3);
    Wig->HSplitPole = HSplitPole;
    Wig->VSplitPole = VSplitPole;
    Wig->zStartH = 0.0;
    Wig->zEndH = Ltot;
    Wig->zStartV = 0.0;
    Wig->zEndV = Ltot;
    /*----------------------------------------------------------*/

    kw = TWOPI/(Wig->Lw);
    Wig->Zw = 0.0;
    Wig->Aw= 1.0e-9*C0/__E0/TWOPI * Lw * Bmax;
    tmppr = By;
    for (i = 0; i < NHharm; i++) {
        tmppr++;
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
        Wig->HCw_raw[i] = 0.0;
        Wig->Hkx[i] = 0.0;
        Wig->Hky[i] = 0.0;
        Wig->Hkz[i] = 0.0;
        Wig->Htz[i] = 0.0;
    }
    for (i = NVharm; i< WHmax; i++) {
        Wig->VCw_raw[i] = 0.0;
        Wig->Vkx[i] = 0.0;
        Wig->Vky[i] = 0.0;
        Wig->Vkz[i] = 0.0;
        Wig->Vtz[i] = 0.0;
    }
}

static double sinc(double x)
{
  double x2, result;
/* Expand sinc(x) = sin(x)/x to x^8 */
  x2 = x*x;
  result = 1e0 - x2/6e0*(1e0 - x2/20e0 *(1e0 - x2/42e0*(1e0-x2/72e0) ) );
  return result;
}

static void GWigAx(struct gwig *pWig, double *Xvec, double *pax, double *paxpy)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sxkx, chx, shx;
  double cy, sy, chy, shy, sz;
  double ax, axpy;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = TWOPI/(pWig->Lw);
  ax   = 0e0;
  axpy = 0e0;

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for (i = 0; i < pWig->NHharm; i++) {
    double HCw = pWig->HCw_raw[i] * pWig->Aw / pWig->Po;
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    cx  = cos(kx*x);
    chy = cosh(ky * y);
    sz  = sin(kz * z + tz);
    ax  = ax + HCw*(kw/kz)*cx*chy*sz;

    shy = sinh(ky * y);
    if (fabs(kx/kw) > GWIG_EPS) {
      sxkx = sin(kx * x)/kx;
    } else {
      sxkx = x*sinc(kx*x);
    }

    axpy = axpy + HCw*(kw/kz)*ky*sxkx*shy*sz;
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
    sz  = sin(kz * z + tz);
    ax  = ax + VCw*(kw/kz)*(ky/kx)*shx*sy*sz;

    chx = cosh(kx * x);
    cy  = cos(ky * y);
    axpy = axpy + VCw*(kw/kz)* pow(ky/kx,2) *chx*cy*sz;
  }

  *pax   = ax;
  *paxpy = axpy;
}

static void GWigAy(struct gwig *pWig, double *Xvec, double *pay, double *paypx)
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sx, chx, shx;
  double cy, syky, chy, shy, sz;
  double ay, aypx;

  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;

  kw   = TWOPI/(pWig->Lw);
  ay   = 0e0;
  aypx = 0e0;

  /* Horizontal Wiggler: note that one potentially could have: kx=0 */
  for (i = 0; i < pWig->NHharm; i++){
    double HCw = pWig->HCw_raw[i] * pWig->Aw / pWig->Po;
    kx = pWig->Hkx[i];
    ky = pWig->Hky[i];
    kz = pWig->Hkz[i];
    tz = pWig->Htz[i];

    sx = sin(kx * x);
    shy = sinh(ky * y);
    sz  = sin(kz * z + tz);
    ay  = ay + HCw*(kw/kz)*(kx/ky)*sx*shy*sz;

    cx  = cos(kx * x);
    chy = cosh(ky * y);

    aypx = aypx + HCw*(kw/kz)*pow(kx/ky,2) * cx*chy*sz;
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
    sz  = sin(kz * z + tz);
    ay  = ay + VCw*(kw/kz)*chx*cy*sz;

    shx = sinh(kx * x);
    if (fabs(ky/kw) > GWIG_EPS) {
      syky  = sin(ky * y)/ky;
    } else {
      syky = y * sinc(ky * y);
    }
    aypx = aypx + VCw*(kw/kz)* kx*shx*syky*sz;
   }

  *pay = ay;
  *paypx = aypx;
}

static void GWigGauge(struct gwig *pWig, double *X, int flag)
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

static void GWigMap_2nd(struct gwig *pWig, double *X, double dl)
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

static void print66(const char *title, double *m)
{
    atPrintf("\n%s\n", title);
    for (int i=0; i<6; i++) {
        int k=6*i;
        atPrintf("%.12g, %.12g, %.12g, %.12g, %.12g, %.12g\n", m[k+0], m[k+1], m[k+2], m[k+3], m[k+4], m[k+5]);
    }
}
