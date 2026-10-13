#ifndef ATVELOCITY_H
#define ATVELOCITY_H

#include <math.h>

/*  Non-relativistic tracking

    The 6th coordinate is ct = beta*c*dt, canonically conjugate to dp = dP/P0,
    where beta is the velocity of the particle. It reduces to the path
    lengthening for ultra-relativistic particles.

    Integrators compute only the geometric path lengthening. In an element of
    length L where dp is constant, the time of flight adds L*(1 - beta/beta0).
    It is applied by the tracking loop in two halves, at the entrance and at
    the exit of each element, which is exact for drift-kick-drift cavities.

    gamma0 <= 0 means ultra-relativistic tracking: no correction.
*/

static void atVelocityShift(double *r_in, int num_particles, double halflength,
                            double gamma0, double betgam0)
/*  Add halflength*(1 - beta/beta0) to ct, computed as
    -dp*(2+dp)/(gamma*(gamma+(1+dp)*gamma0)) to avoid cancellation
*/
{
    int c;
    if ((gamma0 <= 0.0) || (halflength == 0.0)) return;
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) default(none) \
    shared(r_in,num_particles,halflength,gamma0,betgam0) private(c)
    for (c = 0; c < num_particles; c++) {
        double *r6 = r_in+c*6;
        double p = 1.0 + r6[4];
        double betgam = betgam0*p;
        double gamma = sqrt(1.0 + betgam*betgam);
        r6[5] -= halflength*r6[4]*(1.0 + p)/(gamma*(gamma + p*gamma0));
    }
}

#endif /*ATVELOCITY_H*/
