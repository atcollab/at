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

   r6[1] -= L * (ReSum + irho*B[1]*(x*x-0.5*y*y));
   r6[3] += L * (ImSum + irho*B[1]*x*y);
}
