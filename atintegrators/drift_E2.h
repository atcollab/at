#ifdef DIFFUSION
#error "drift_E2 does not compute the diffusion matrix"
#endif

/* the pseudo-drift element described by Hamiltonian H1 = (1+hx) (px^2+py^2)/2(1+delta),     */
static void drift(double* r6, double L, double h, double *bdiff)
{
    double p_norm = 1.0 / (1.0+r6[4]);
    double px = r6[1];
    double py = r6[3];
	double hs = h*L;
	double x=r6[0];

	r6[0] += (1.0+h*x)*L*p_norm*px + 1.0/4.0*hs*L*(px*px-py*py)*p_norm*p_norm;
	r6[1] -= hs*(px*px+py*py)*p_norm/2.0;

	r6[2] += (1.0+h*x)*L*p_norm*py*(1.0+px*hs/2.0);
	r6[5] += (1.0+h*x)*L*p_norm*p_norm/2.0*(px*px+py*py);
}
