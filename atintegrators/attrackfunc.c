#include <math.h>
#include "atelem.c"
#include "driftkickrad.c"

static void relativistic_kick(double *r6, double nv, double freq, double lag, double philag,
                              double tau_ref)
{
    double tau_rel = (r6[5]-lag)/C0;
    r6[4] += -nv*sin(TWOPI*freq*(tau_rel + tau_ref) - philag);
}

static void non_relativistic_kick(double *r6, double nv, double freq, double lag, double philag,
                                  double tau_ref, double gamma0, double betgam0)
/*  Exact energy gain: with E/P0c = sqrt((1+dp)^2 + 1/betgam0^2), the energy
    gain normalised to P0c is nv/beta0. ct = beta*c*dt is rescaled to keep
    the time continuous.
*/
{
    double a2 = 1.0/(betgam0*betgam0);
    double beta0 = betgam0/gamma0;
    double pi = 1.0 + r6[4];
    double ei = sqrt(pi*pi + a2);
    double betai = pi/ei;
    double tau_rel = (r6[5]/betai - lag/beta0)/C0;
    double de = -nv/beta0*sin(TWOPI*freq*(tau_rel + tau_ref) - philag);
    double ef = ei + de;
    double pf = sqrt(ef*ef - a2);
    double betaf = pf/ef;
    /*  pf-pi computed without cancellation */
    r6[4] += de*(ef + ei)/(pf + pi);
    /*  ensure continuity of dt = ct / beta / clight */
    r6[5] *= betaf/betai;
}

void trackRFCavity(double *r_in, double le, double nv, double freq, double h, double lag, double philag,
                  int nturn, double T0, double gamma0, int num_particles)
/* le - physical length
   nv - peak voltage (V) normalized to the design enegy (eV)
   gamma0 - Lorentz factor of the reference particle, <= 0 for ultra-relativistic tracking
   r is a 6-by-N matrix of initial conditions reshaped into
   1-d array of 6*N elements
*/
{
    int c;
    double tau_ref = -(h/freq-T0)*nturn;

    /* If nv is 0 and length is 0, then skip this whole loop (good for passive rf cavities
        anyway if there is a cavity length, we have to loop through the particles
    */

    if (le == 0) {
        if (nv != 0) {
            for (c = 0; c<num_particles; c++) {
                double *r6 = r_in+c*6;
                if(!atIsNaN(r6[0])) {
                    if (gamma0 > 0.0) {
                        double betgam0 = sqrt(gamma0*gamma0 - 1.0);
                        non_relativistic_kick(r6, nv, freq, lag, philag, tau_ref, gamma0, betgam0);
                    }
                    else
                        relativistic_kick(r6, nv, freq, lag, philag, tau_ref);
                }
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
                if(nv!=0.0) {
                    if (gamma0 > 0.0) {
                        double betgam0 = sqrt(gamma0*gamma0 - 1.0);
                        non_relativistic_kick(r6, nv, freq, lag, philag, tau_ref, gamma0, betgam0);
                    }
                    else
                        relativistic_kick(r6, nv, freq, lag, philag, tau_ref);
                }
                /* Propagate through a drift equal to half cavity length */
                drift6(r6, halflength);
            }
        }
    }
}
