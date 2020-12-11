/*
 *----------------------------------------------------------------------------
 * Modification Log:
 * -----------------
 * .04  2003-04-29      YK Wu, Duke University, wu@fel.duke.edu 
 *			using scientific notation for constants. 
 *                      Checked with TRACY pascal code.
 *                      Computing differential pathlength only.
 *
 * .03  2003-04-28      YK Wu, Duke University, wu@fel.duke.edu 
 *			Convert to C code and cross-checked with the pascal version;
 *
 * .02  2001-12-xx      Y. K. Wu, Duke University, wu@fel.duke.edu
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

void GWigGauge(struct gwigR *pWig, double *X, int flag);
void GWigPass_2nd(struct gwigR *pWig, double *X);
void GWigPass_4th(struct gwigR *pWig, double *X);
void GWigMap_2nd(struct gwigR *pWig, double *X, double dl);
void GWigAx(struct gwigR *pWig, double *Xvec, double *pax, double *paxpy);
void GWigAy(struct gwigR *pWig, double *Xvec, double *pay, double *paypx);
void GWigRadiationKicks(struct gwigR *pWig, double *X, double *Bxyz, double dl);
void GWigB(struct gwigR *pWig, double *Xvec, double *B);
double sinc(double x );

/* This function appears to be unused. */
void GWigGauge(struct gwigR *pWig, double *X, int flag)
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


void GWigPass_2nd(struct gwigR *pWig, double *X) 
{
  int    i, Nstep;
  double dl;
  double B[2];
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


void GWigPass_4th(struct gwigR *pWig, double *X)
{

  const double x1 = 1.3512071919596576340476878089715e0;
  const double x0 =-1.7024143839193152680953756179429e0;

  int    i, Nstep;
  double dl, dl1, dl0;
  double B[2];
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
  for (i = 1; i <= Nstep; i++ ) {
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


void GWigMap_2nd(struct gwigR *pWig, double *X, double dl) 
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


void GWigAx(struct gwigR *pWig, double *Xvec, double *pax, double *paxpy) 
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
  for (i = 0; i < pWig->NVharm; i++ ) {
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


void GWigAy(struct gwigR *pWig, double *Xvec, double *pay, double *paypx)
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
  for ( i = 0; i < pWig->NHharm; i++ ){
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


double sinc(double x)
{
  double x2, result;
/* Expand sinc(x) = sin(x)/x to x^8 */
  x2 = x*x;
  result = 1e0 - x2/6e0*(1e0 - x2/20e0 *(1e0 - x2/42e0*(1e0-x2/72e0) ) );
  return result;
}



void GWigB(struct gwigR *pWig, double *Xvec, double *B) 
/* Compute magnetic field at particle location.
 * Added by M. Borland, August 2007.
 */
{
  int    i;
  double x, y, z;
  double kx, ky, kz, tz, kw;
  double cx, sx, chx, shx;
  double cy, sy, chy, shy;
  double cz;
  /* B0 is a reserved symbol on MacOS, defined in termios.h */
  double _B0;
  
  x = Xvec[0];
  y = Xvec[2];
  z = pWig->Zw;
  
  kw   = 2e0*PI/(pWig->Lw);

  B[0] = 0;
  B[1] = 0;

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
        
        /* Accumulate field values in user-supplied array (Bx, By) */
        B[0] += _B0*pWig->HCw_raw[i]*kx/ky*sx*shy*cz;
        B[1] -= _B0*pWig->HCw_raw[i]*cx*chy*cz;
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
        
        B[0] -= _B0*pWig->HCw_raw[i]*kx/ky*shx*sy*cz;
        B[1] -= _B0*pWig->HCw_raw[i]*chx*cy*cz;
      }
    }
  }
  
  if (pWig->NVharm && z>=pWig->zStartV && z<=pWig->zEndV) {
    _B0 = pWig->PB0;
    if (!pWig->VSplitPole) {
      /* Normal Vertical Wiggler: note that one potentially could have: ky=0 */
      for (i = 0; i < pWig->NVharm; i++ ) {
        kx = pWig->Vkx[i];
        ky = pWig->Vky[i];
        kz = pWig->Vkz[i];
        tz = pWig->Vtz[i];

        shx = sinh(kx * x);
        chx = cosh(kx * x);
        sy  = sin(ky * y);
        cy  = cos(ky * y);
        cz  = cos(kz*z + tz);

        /* Accumulate field values in user-supplied array (Bx, By) */
        B[0] += _B0*pWig->VCw_raw[i]*chx*cy*cz;
        B[1] -= _B0*pWig->VCw_raw[i]*ky/kx*shx*sy*cz;
      }
    } else {
      /* Split-pole Vertical Wiggler: note that one potentially could have: kx=0 (caught in main routine) */
      for (i = 0; i < pWig->NVharm; i++ ) {
        kx = pWig->Vkx[i];
        ky = pWig->Vky[i];
        kz = pWig->Vkz[i];
        tz = pWig->Vtz[i];

        sx  = sin(kx * x);
        cx  = cos(kx * x);
        shy = sinh(ky * y);
        chy = cosh(ky * y);
        cz  = cos(kz*z + tz);

        /* Accumulate field values in user-supplied array (Bx, By) */
        B[0] += _B0*pWig->VCw_raw[i]*cx*chy*cz;
        B[1] += _B0*pWig->VCw_raw[i]*ky/kx*sx*shy*cz;
      }
    }
  }
}


void GWigRadiationKicks(struct gwigR *pWig, double *X, double *Bxy, double dl)
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

