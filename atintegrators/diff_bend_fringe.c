
static void diff_bend_fringe(double *r6, double inv_rho, double edge_angle, double fint,
                             double gap, int method, double sign, double *bdiff) {
  /*     method 0 no fringe field
   *     method 1 legacy version Brown First Order
   *     method 2 SOLEIL close to second order of Brown
   *     method 3 THOMX
   */
  double psi, fx, fy;
  double p_norm = 1.0 / (1.0+r6[4]);
  /* Fringe field correction */
  if ((fint == 0.0) || (gap == 0.0) || (method == 0))
    psi = 0.0;
  else {
    register double sedge = sin(edge_angle);
    register double cedge = cos(edge_angle);
    psi = inv_rho * gap * fint * (1.0 + sedge*sedge) / cedge;
  }

  /* Edge angle focusing */
  fx = inv_rho * tan(edge_angle);
  if (method == 1)
    fy = inv_rho * tan(edge_angle - psi*p_norm);
  else if (method == 2)
    fy = inv_rho * tan(edge_angle - psi*p_norm) * p_norm;
  else if (method == 3)
    fy = inv_rho * tan(edge_angle - psi + sign*r6[1]*p_norm);
  else /* fall back to legacy version */
    fy = inv_rho * tan(edge_angle - psi*p_norm);

  /*  Propagate B */
  if (bdiff) {
      for (int m = 0; m < 6; m++) {
        bdiff[1 + 6*m] += fx * bdiff[6*m];
        bdiff[3 + 6*m] -= fy * bdiff[2 + 6*m];
      }
      if (fint > 0 && gap > 0)
        for (int m = 0; m < 6; m++)
          bdiff[3 + 6*m] -= bdiff[4 + 6*m] * r6[2] * (inv_rho*inv_rho + fy*fy) * psi * p_norm * p_norm / inv_rho;

      for (int m = 0; m < 6; m++) {
        bdiff[m + 6*1] += fx * bdiff[m + 6*0];
        bdiff[m + 6*3] -= fy * bdiff[m + 6*2];
      }
      if (fint > 0 && gap > 0)
        for (int m = 0; m < 6; m++)
          bdiff[m + 6*3] -= bdiff[m + 6*4] * r6[2] * (inv_rho*inv_rho + fy*fy) * psi * p_norm * p_norm / inv_rho;
  }
  /*  Propagate particle */

  r6[1] += r6[0] * fx;
  r6[3] -= r6[2] * fy;
}
