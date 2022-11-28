#include <math.h>
#include "atelem.c"
#include "driftkickrad.c"

#define TWOPI  6.28318530717959
#define C0  	2.99792458e8 

void trackRFCavity(double *r_in, double le, double nv, double freq, double h, double lag, double philag,
                  int nturn, double T0, int num_particles)
/* le - physical length
   nv - peak voltage (V) normalized to the design enegy (eV)
   r is a 6-by-N matrix of initial conditions reshaped into
   1-d array of 6*N elements
*/
{
    int c;

    if (le == 0) {
        for (c = 0; c<num_particles; c++) {
            double *r6 = r_in+c*6;
            if(!atIsNaN(r6[0]))
                r6[4] += -nv*sin(TWOPI*freq*((r6[5]-lag)/C0 - (h/freq-T0)*nturn) - philag);
        }
    }
    else {
        double halflength = le/2;
        for (c = 0;c<num_particles;c++) {
            double *r6 = r_in+c*6;
            if(!atIsNaN(r6[0]))  {
                /* Propagate through a drift equal to half cavity length */
                drift6(r6, halflength);
                /* Longitudinal momentum kick */
                r6[4] += -nv*sin(TWOPI*freq*((r6[5]-lag)/C0 - (h/freq-T0)*nturn) - philag);
                /* Propagate through a drift equal to half cavity length */
                drift6(r6, halflength);
            }
        }
    }
}
