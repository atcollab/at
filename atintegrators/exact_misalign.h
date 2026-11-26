#ifndef EXACT_MISALIGN_H
#define EXACT_MISALIGN_H

/*
 * Exact misalignment transformations at the faces of an element.
 *
 * M. Venturini and T. Hellert, "Particle-ray tracing through misaligned
 * lattice elements in accelerators", Nucl. Instrum. Methods A,
 * https://doi.org/10.1016/j.nima.2022.167487
 *
 * at.lattice.transform_elem builds the linearised face maps of that paper
 * (Eqs. 8-11 and 27-29) and stores them in T1, R1 (r <- R1 (r + T1)) and
 * R2, T2 (r <- R2 r + T2). With ExactMisalign=1, the rotation R and the
 * displacement d of each face are recovered from these matrices, and the
 * exact map of Secs. 3.1, 4.1 and Appendix A is applied instead:
 *
 *      (X, Y, Z)    = R^T (x, y, 0) - u,       u = R^T d
 *      (Px, Py, Pz) = R^T (px, py, pz),        pz = sqrt((1+delta)^2 - px^2 - py^2)
 *
 * followed by an exact drift over L = -Z:
 *
 *      x = X + L Px/Pz,   y = Y + L Py/Pz,   ct += L (1+delta)/Pz
 *
 * Recovering the geometry: the px and py rows of the linear map, in the
 * (px, py, delta) columns, are the first two columns of R. Its constant term
 * c (R1 T1 at the entrance, T2 at the exit) is the exact map of r = 0:
 *
 *      c_ct = u_z / R_zz,   c_x = u_z R_zx / R_zz - u_x,   c_y = u_z R_zy / R_zz - u_y
 *
 * The frame tilt (tilt_frame) is a rotation about the Z axis of the face,
 * already included in R1 and R2, and so in R.
 *
 * The geometry is stored as 12 doubles: R row by row, then u.
 */

/* Face geometry from the linear map M (column-major 6x6, NULL for identity)
   and its constant term c (NULL for zero). Returns 0 if M is not a rigid
   misalignment matrix. */
static int exact_misalign_geometry(const double *M, const double *c, double *g)
{
    double *R = g;
    double *u = g + 9;
    double c0 = c ? c[x_] : 0.0;
    double c2 = c ? c[y_] : 0.0;
    double c5 = c ? c[ct_] : 0.0;

    if (M) {
        R[0] = M[px_ + 6*px_]; R[3] = M[px_ + 6*py_]; R[6] = M[px_ + 6*delta_];
        R[1] = M[py_ + 6*px_]; R[4] = M[py_ + 6*py_]; R[7] = M[py_ + 6*delta_];
    }
    else {
        R[0] = 1.0; R[3] = 0.0; R[6] = 0.0;
        R[1] = 0.0; R[4] = 1.0; R[7] = 0.0;
    }
    /* Third column: cross product of the first two */
    R[2] = R[3]*R[7] - R[6]*R[4];
    R[5] = R[6]*R[1] - R[0]*R[7];
    R[8] = R[0]*R[4] - R[3]*R[1];

    /* The two columns must be orthonormal */
    double n0 = R[0]*R[0] + R[3]*R[3] + R[6]*R[6];
    double n1 = R[1]*R[1] + R[4]*R[4] + R[7]*R[7];
    double p01 = R[0]*R[1] + R[3]*R[4] + R[6]*R[7];
    if (fabs(n0 - 1.0) > 1.0e-10 || fabs(n1 - 1.0) > 1.0e-10 || fabs(p01) > 1.0e-10)
        return 0;

    u[2] = c5 * R[8];
    u[0] = u[2] * R[6] / R[8] - c0;
    u[1] = u[2] * R[7] / R[8] - c2;
    return 1;
}

/* Geometry of both faces in geometry[24]. E1 (E2) points to the entrance
   (exit) part, or is NULL if the face is not transformed. Returns 0 if R1 or
   R2 is not a rigid misalignment matrix. */
static int exact_misalign_init(const double *R1, const double *T1,
        const double *R2, const double *T2,
        double *geometry, double **E1, double **E2)
{
    *E1 = NULL;
    *E2 = NULL;
    if (R1 || T1) {
        double c[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        if (T1) {
            for (int i = 0; i < 6; i++) c[i] = T1[i];
            if (R1) ATmultmv(c, R1);
        }
        if (!exact_misalign_geometry(R1, c, geometry)) return 0;
        *E1 = geometry;
    }
    if (R2 || T2) {
        if (!exact_misalign_geometry(R2, T2, geometry + 12)) return 0;
        *E2 = geometry + 12;
    }
    return 1;
}

#ifdef DIFFUSION
/* Jacobian of exact_misalign at r6, column-major as R1, R2 */
static void exact_misalign_jacobian(const double *r6, const double *g, double *M)
{
    const double *R = g;
    const double *u = g + 9;
    double dp1 = 1.0 + r6[delta_];
    double pz = sqrt(dp1*dp1 - r6[px_]*r6[px_] - r6[py_]*r6[py_]);
    double Px = R[0]*r6[px_] + R[3]*r6[py_] + R[6]*pz;
    double Py = R[1]*r6[px_] + R[4]*r6[py_] + R[7]*pz;
    double Pz = R[2]*r6[px_] + R[5]*r6[py_] + R[8]*pz;
    double L = u[2] - R[2]*r6[x_] - R[5]*r6[y_];
    /* Derivatives with respect to (x, px, y, py, delta, ct) */
    double gpz[6] = {0.0, -r6[px_]/pz, 0.0, -r6[py_]/pz, dp1/pz, 0.0};
    double gdp1[6] = {0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
    double gX[6] = {R[0], 0.0, R[3], 0.0, 0.0, 0.0};
    double gY[6] = {R[1], 0.0, R[4], 0.0, 0.0, 0.0};
    double gL[6] = {-R[2], 0.0, -R[5], 0.0, 0.0, 0.0};
    double gPx[6], gPy[6], gPz[6];

    for (int j = 0; j < 6; j++) {
        gPx[j] = R[6]*gpz[j];
        gPy[j] = R[7]*gpz[j];
        gPz[j] = R[8]*gpz[j];
    }
    gPx[px_] += R[0]; gPx[py_] += R[3];
    gPy[px_] += R[1]; gPy[py_] += R[4];
    gPz[px_] += R[2]; gPz[py_] += R[5];

    for (int j = 0; j < 6; j++) {
        M[x_ + 6*j] = gX[j] + gL[j]*Px/Pz + L*(gPx[j]*Pz - Px*gPz[j])/(Pz*Pz);
        M[px_ + 6*j] = gPx[j];
        M[y_ + 6*j] = gY[j] + gL[j]*Py/Pz + L*(gPy[j]*Pz - Py*gPz[j])/(Pz*Pz);
        M[py_ + 6*j] = gPy[j];
        M[delta_ + 6*j] = gdp1[j];
        M[ct_ + 6*j] = gL[j]*dp1/Pz + L*(gdp1[j]*Pz - dp1*gPz[j])/(Pz*Pz);
    }
    M[ct_ + 6*ct_] += 1.0;
}
#endif /*DIFFUSION*/

static void exact_misalign(double *r6, const double *g)
{
    const double *R = g;
    const double *u = g + 9;
    double dp1 = 1.0 + r6[delta_];
    double pz = sqrt(dp1*dp1 - r6[px_]*r6[px_] - r6[py_]*r6[py_]);
    double X = R[0]*r6[x_] + R[3]*r6[y_] - u[0];
    double Y = R[1]*r6[x_] + R[4]*r6[y_] - u[1];
    double L = u[2] - R[2]*r6[x_] - R[5]*r6[y_];
    double Px = R[0]*r6[px_] + R[3]*r6[py_] + R[6]*pz;
    double Py = R[1]*r6[px_] + R[4]*r6[py_] + R[7]*pz;
    double Pz = R[2]*r6[px_] + R[5]*r6[py_] + R[8]*pz;

    r6[x_] = X + L*Px/Pz;
    r6[px_] = Px;
    r6[y_] = Y + L*Py/Pz;
    r6[py_] = Py;
    r6[ct_] += L*dp1/Pz;
}

#endif /*EXACT_MISALIGN_H*/
