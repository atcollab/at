
static void drift(double* r, double NormL)

/*   NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
 in the loop if momentum deviation (delta) does not change
 such as in 4-th order symplectic integrator w/o radiation
 */

{
   r[0] += NormL*r[1];
   r[2] += NormL*r[3];
   r[5] += NormL*(r[1]*r[1]+r[3]*r[3])/(2*(1+r[4]));
}
