static void kick(double* r6, double A0, double B0, const double* A, const double* B, int max_order, double L, double irho)
{
   double ReSum = B[max_order];
   double ImSum = A[max_order];
   double ReSumTemp;
   double x = r6[0];
   double y = r6[2];

   /* recursively calculate the local transverse magnetic field */
   for (int i=max_order-1; i>=0; i--) {
     ReSumTemp = ReSum*x - ImSum*y + B[i];
     ImSum = ImSum*x +  ReSum*y + A[i];
     ReSum = ReSumTemp;
   }
   ReSum += B0;
   ImSum += A0;

   r6[1] -= L * ReSum;
   r6[3] += L * ImSum;

   /* Curvature corrections.  In the curved frame the multipole potential
      carries the metric factor (1 + irho*x), which produces one correction per
      field order.  The design dipole is handled by the exact bend propagator,
      so only the *error* dipole B[0] (plus the corrector kick B0) contributes
      here: H = 1/2 * irho * b0 * x^2  (MAD8 physics manual eq. 5.15). */
   r6[1] -= L * irho * (B[0] + B0) * x;

   if (max_order >= 1) {  /* K1h correction */
     r6[1] -= L * irho * B[1] * (x*x-0.5*y*y);
     r6[3] += L * irho * B[1] * x * y;
   }
}
