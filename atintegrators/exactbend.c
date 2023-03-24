
static void exact_bend(double *r6, double irho, double L)
{
    double dp1 = 1.0 + r6[delta_];
    double pz = pxyz(dp1, r6[px_], r6[py_]);
    double cs = cos(irho*L);
    double sn = sin(irho*L);
    double px = r6[px_]*cs + (pz-(1.0+r6[x_]*irho))*sn;
    double d2 = pxyz(dp1, 0.0, r6[py_]);
    double dasin = L + (asin(r6[px_]/d2) - asin(px/d2))/irho;
    double x = (pxyz(dp1,px,r6[py_]) - (pz-(1.0+r6[x_]*irho))*cs + r6[px_]*sn - 1.0)/irho;
    double dy = r6[py_]*dasin;
    double dct = dp1*dasin;
    
    r6[x_] = x;
    r6[px_] = px;
    r6[y_] += dy;
    r6[ct_] += dct;
}