#include <math.h>
#include "atelem.c"
#include "driftkickrad.c"

void trackRFCavity(double *r_in, double le, double nv, double freq, double h, double lag, double philag,
                  int nturn, double T0, int num_particles)
/* le - physical length
   nv - peak voltage (V) normalized to the design enegy (eV)
   r is a 6-by-N matrix of initial conditions reshaped into
   1-d array of 6*N elements
*/
{
    int c;

    /* If nv is 0 and length is 0, then skip this whole loop (good for passive rf cavities
        anyway if there is a cavity length, we have to loop through the particles
    */

    if (le == 0) {
        if (nv != 0) {
            for (c = 0; c<num_particles; c++) {
                double *r6 = r_in+c*6;
                if(!atIsNaN(r6[0]))
                    r6[4] += -nv*sin(TWOPI*freq*((r6[5]-lag)/C0 - (h/freq-T0)*nturn) - philag);
            }
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
                if(nv!=0.0) r6[4] += -nv*sin(TWOPI*freq*((r6[5]-lag)/C0 - (h/freq-T0)*nturn) - philag);
                /* Propagate through a drift equal to half cavity length */
                drift6(r6, halflength);
            }
        }
    }
}
