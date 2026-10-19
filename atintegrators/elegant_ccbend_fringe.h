/************************************************************************* \
 * Copyright (c) 2018 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2018 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution.
\*************************************************************************/

/* Modified for AT, 2026: isolated the two end-map functions, made them static;
 * formulas and USE_NEW_SWITCH_CODE=0 retained verbatim.
 * License agreement: elegant.LICENSE (distributed with this source).
 * Original ccbend.c SHA256: eb4f0dd128ff74efdb880cedf1eed9f11ac5283375817a4d9161ce5060147039
 */
#ifndef AT_ELEGANT_CCBEND_FRINGE_H
#define AT_ELEGANT_CCBEND_FRINGE_H
#define sqr(x) ((x)*(x))
#define USE_NEW_SWITCH_CODE 0
static void switchRbendPlane(double **particle, long n_part, double alpha, double po)
/* transforms the reference plane to one that is at an angle alpha relative to the
 * initial plane. use alpha=theta/2, where theta is the total bend angle.
 */
{
#if USE_NEW_SWITCH_CODE
  long i;
  double cos_alpha, tan_alpha, magnet_s;

  tan_alpha = tan(alpha);
  cos_alpha = cos(alpha);

  for (i = 0; i < n_part; i++) {
    magnet_s = particle[i][0] * tan_alpha / (1.0 - particle[i][1] * tan_alpha);
    particle[i][0] = (particle[i][0] + particle[i][1] * magnet_s) / cos_alpha;
    particle[i][2] = particle[i][2] + particle[i][3] * magnet_s;

    particle[i][4] += sqrt(1.0 + sqr(particle[i][1]) + sqr(particle[i][3])) * magnet_s;

    particle[i][3] = particle[i][3] / (cos_alpha * (1.0 - particle[i][1] * tan_alpha));
    particle[i][1] = (particle[i][1] + tan_alpha) / (1.0 - particle[i][1] * tan_alpha);
  }
#else
  long i;
  double s, *coord, sin_alpha, cos_alpha, tan_alpha;
  double d, qx0, qy0, qz0, qx, qy, qz;

  tan_alpha = tan(alpha);
  cos_alpha = cos(alpha);
  sin_alpha = sin(alpha);

  for (i = 0; i < n_part; i++) {
    coord = particle[i];
    s = coord[0] * tan_alpha / (1 - coord[1] * tan_alpha);
    coord[0] = (coord[0] + coord[1] * s) / cos_alpha;
    coord[2] = (coord[2] + coord[3] * s);
    coord[4] += s;
    d = sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
    qx0 = coord[1] * (1 + coord[5]) / d;
    qy0 = coord[3] * (1 + coord[5]) / d;
    qz0 = (1 + coord[5]) / d;
    qx = qx0 * cos_alpha + qz0 * sin_alpha;
    qy = qy0;
    qz = -qx0 * sin_alpha + qz0 * cos_alpha;
    coord[1] = qx / qz;
    coord[3] = qy / qz;
  }
#endif
}

static void verticalRbendFringe(
                         double **particle, long n_part,
                         double alpha,        // edge angle relative to beam path
                         double rho0,         // bending radius (always positive)
                         double K1,           // interior gradient
                         double K2,           // interior sextupole
                         double gK,           // conventional soft-edge parameter
                         double *fringeIntKn, // fringe integrals: K0, I0, K2, I1, K4, K5, K6, K7
                         short angleSign,     // -1 or 1
                         short isExit,
                         short fringeModel,
                         short order) {

  if (fringeModel == -1)
    return;

  if (fringeModel == 0) {
    // old method
    long i;
    double c, d, e;
    if (order < 1)
      return;
    c = d = e = 0;
    if (gK != 0)
      alpha -= gK / fabs(rho0) / cos(alpha) * (1 + sqr(sin(alpha)));
    c = sin(alpha) / rho0;
    if (order > 1)
      d = sin(alpha) * K1;
    if (order > 2)
      e = sin(alpha) * K2 / 2;
    for (i = 0; i < n_part; i++) {
      double x, y;
      x = particle[i][0];
      y = particle[i][2];
      particle[i][3] -= y * (c + d * x + e * (x * x - y * y)) / (1 + particle[i][5]);
    }
  } else {
    double x1, px1, y1, py1, tau1, delta;
    double x2, px2, y2, py2, tau2;

    double intK0, intK2, intK4, intK5, intK6, intK7, intI0, intI1;
    double invRhoPlus, invRhoMinus, K1plus, K1minus;
    double tant, sect, sect3, sint, temp;
    double focX0, /* focXd, */ focY0, focYd, invP;
    double dispX, kickPx, expT;

    long i;
    if (isExit)
      alpha = -alpha;

    // Per R. Lindberg, some of the fringe integrals change sign with the sign of the bending angle.
    intK0 = fringeIntKn[0] * angleSign;
    intK2 = fringeIntKn[2];
    intK4 = fringeIntKn[4] * angleSign;
    intK5 = fringeIntKn[5] * angleSign;
    intK6 = fringeIntKn[6] * angleSign;
    intK7 = fringeIntKn[7];
    // shamelessly packed the In integrals in the slots unused by the Kn's
    intI0 = fringeIntKn[1];
    intI1 = fringeIntKn[3];

    if (isExit) {
      // exit fringe
      K1minus = K1;
      K1plus = 0;
      invRhoMinus = 1 / rho0;
      invRhoPlus = 0;
    } else {
      // entrance fringe
      K1plus = K1;
      K1minus = 0;
      invRhoPlus = 1 / rho0;
      invRhoMinus = 0;
    }

    tant = tan(alpha);
    sint = sin(alpha);
    sect = 1.0 / cos(alpha);
    sect3 = sect * sect * sect;
    focX0 = -tant * intK5 - 0.5 * intI0 * (2.0 - tant * tant);
    /* focXd = sect3 * intK7; */
    focY0 = -tant * (invRhoPlus - invRhoMinus) + tant * sect * sect * intK5 + 0.5 * intI0 * (2.0 + tant * tant);
    focYd = sect3 * ((1.0 + sint * sint) * intK2 - intK7);
    for (i = 0; i < n_part; i++) {
      x1 = particle[i][0];
      temp = sqrt(1.0 + sqr(particle[i][1]) + sqr(particle[i][3]));
      px1 = (1.0 + particle[i][5]) * particle[i][1] / temp;
      y1 = particle[i][2];
      py1 = (1.0 + particle[i][5]) * particle[i][3] / temp;
      tau1 = -particle[i][4] + sint * x1;
      delta = particle[i][5];
      px1 = px1 - (1.0 + delta) * sint;
      invP = 1.0 / (1.0 + delta);

      /** Symplectic update for Linear and Zeroth Order terms **/
      temp = sect * intI1 * invP;
      if (fabs(temp) < 1.0e-20)
        temp = 1.0e-20;
      dispX = -sect * intK0 * invP;
      kickPx = -tant * (intI1 + 0.5 * intK4);
      expT = exp(temp);

      x2 = expT * x1 + dispX * (expT - 1.0) / temp;
      y2 = y1 / expT;
      px2 = px1 / expT + kickPx * (1.0 - 1.0 / expT) / temp + focX0 * (x1 * (expT - 1.0 / expT) / (2.0 * temp) + dispX * (expT + 1.0 / expT - 2.0) / (2.0 * temp * temp));
      py2 = expT * py1 + y1 * (focY0 + focYd * invP) * (expT - 1.0 / expT) / (2.0 * temp);
      tau2 = tau1 + invP * (kickPx * (temp * x1 + dispX) * (1.0 + temp - expT) / (temp * temp) - px1 * dispX - (px1 * x1 - py1 * y1) * temp + (0.5 * focYd * invP + focY0 * (1.0 / (expT * expT) - 1.0 + 2.0 * temp) / (4.0 * temp)) * y1 * y1 - x1 * x1 * focX0 * (expT * expT - 1.0 - 2.0 * temp) / (4.0 * temp) - x1 * dispX * focX0 * (expT * expT - 2.0 * expT + 1.0) / (2.0 * temp * temp) - dispX * dispX * focX0 * (expT * expT - 4.0 * expT + 2.0 * temp + 3.0) / (4.0 * temp * temp * temp));

      /* x2 += 0.5*sect3*( intK5*(x1*x1 - y1*y1) + y1*y1*(invRhoPlus-invRhoMinus) )*invP;
         y2 += -sect*intK5*x1*y1*invP;
         px2 += -0.25*tant*(K1plus - K1minus)*(x1*x1 + y1*y1) - 0.5*intK6*(x1*x1 - sect*sect*y1*y1)
         + intK5*(sect*py1*y1 - sect3*px1*x1)*invP; // hereit
         py2 += (sect*sect*intK6 - 0.5*tant*(K1plus - K1minus))*x1*y1
         + ( intK5*(sect*py1*x1 + sect3*px1*y1) - sect3*(invRhoPlus - invRhoMinus)*px1*y1 )*invP;
         tau2 += ( sect*intK5*py1*x1*y1
         + sect3*(intK5*px1*(y1*y1-x1*x1) - px1*y1*y1*(invRhoPlus-invRhoMinus)) )*invP*invP; */

      if (order >= 2) {
        /** Symplectic update for second order terms **/
        x1 = x2;
        y1 = y2 * exp(-sect * intK5 * x2 * invP);
        px1 = px2 - (0.5 * intK6 + 0.25 * tant * (K1plus - K1minus)) * x2 * x2 + sect * intK5 * invP * py2 * y2;
        py1 = py2 * exp(sect * intK5 * x2 * invP);
        tau1 = tau2 + sect * intK5 * invP * invP * py2 * x2 * y2;

        x2 = x1 - 0.5 * sect3 * (intK5 - (invRhoPlus - invRhoMinus)) * invP * y1 * y1;
        y2 = y1;
        px2 = px1 + (0.5 * sect * sect * intK6 - 0.25 * tant * (K1plus - K1minus)) * y1 * y1;
        py2 = py1 + sect3 * (intK5 - (invRhoPlus - invRhoMinus)) * invP * px1 * y1 + (sect * sect * intK6 - 0.5 * tant * (K1plus - K1minus)) * x1 * y1;
        tau2 = tau1 + 0.5 * sect3 * (intK5 - (invRhoPlus - invRhoMinus)) * invP * invP * y1 * y1 * (px1 + (0.25 * sect * sect * intK6 - 0.125 * tant * (K1plus - K1minus)) * y1 * y1);

        temp = -0.5 * sect3 * intK5 * invP * x2;
        tau2 += temp * invP * px2 * x2;
        temp = 1.0 + temp;
        x2 = x2 / temp;
        px2 = px2 * temp * temp;
      }

      particle[i][0] = x2;
      particle[i][2] = y2;
      px2 = px2 + (1.0 + delta) * sint;
      temp = sqrt(sqr(1 + delta) - sqr(px2) - sqr(py2));
      particle[i][1] = px2 / temp;
      particle[i][3] = py2 / temp;
      particle[i][4] = -tau2 + sint * x2;
    }
  }
}


#undef sqr
#undef USE_NEW_SWITCH_CODE
#endif
