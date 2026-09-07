#ifndef EXACT_MISALIGN_H
#define EXACT_MISALIGN_H

/*
 * Exact (nonlinear) rigid-body misalignment for the Exact* pass methods.
 *
 * AT normally applies misalignments through the linearised 6x6 R1/R2 matrices
 * and T1/T2 vectors built by at.lattice.transformation.transform_elem().  That
 * linearisation is inconsistent with an "exact" integrator: the geometry it
 * encodes is exact, but its action on the particle is truncated to first order
 * in the transverse momenta.  This header applies the *same* geometry exactly,
 * reproducing Xsuite's track_misalignments.h (xtrack >= 0.109).
 *
 * Conventions
 * -----------
 * Xsuite parameters are obtained from the AT element attributes as
 *
 *      theta (rot_y_rad)          =  yaw
 *      phi   (rot_x_rad)          = -pitch
 *      psi   (rot_s_rad_no_frame) =  tilt
 *      dx, dy, ds                 =  dx, dy, dz
 *
 * (the pitch sign flip is the one already used by at.load.xsuite).
 *
 * The anchor is the distance from the element entrance to the point the
 * transformation is referenced to:  Length/2 for ReferencePoint.CENTRE (the AT
 * default) and 0 for ReferencePoint.ENTRANCE.
 *
 * Longitudinal coordinate: AT stores the *absolute path length* in r6[ct_]
 * inside the Exact* pass methods, whereas Xsuite advances zeta = s - beta0*c*t.
 * The two run in opposite directions, so every longitudinal increment below is
 * the negative of its Xsuite counterpart.  Consequently the two rotations map
 * onto Xsuite's macros with a reversed angle:
 *
 *      Xsuite Y_ROTATE(a)  ==  mis_yrot(r6, -a)
 *      Xsuite X_ROTATE(a)  ==  mis_xrot(r6, -a)
 *      Xsuite S_ROTATE(a)  ==  mis_srot(r6,  a)
 *
 * mis_yrot() is deliberately identical in form to the Yrot() of
 * exactbendfringe.c (Forest 10.26), which fixes the AT sign convention.
 */

#include <math.h>

#ifndef PXYZ
#define PXYZ
static double pxyz(double dp1, double px, double py)
{
  return sqrt(dp1*dp1 - px*px - py*py);
}
#endif /*PXYZ*/

/* ------------------------------------------------------------------ */
/* Elementary exact transformations, in the AT coordinate convention   */
/* ------------------------------------------------------------------ */

static void mis_xy_shift(double *r6, double dx, double dy)
{
    r6[x_] -= dx;
    r6[y_] -= dy;
}

static void mis_s_shift(double *r6, double ds)
{
    /* Longitudinal displacement of the element: the particle covers the real
       path over ds, but the *design* orbit does not advance, so the whole
       trajectory length is path lengthening.  This is the AT image of Xsuite's
           Drift_single_particle_exact(ds); zeta -= ds; s -= ds;
       whose net effect on zeta is -ds*(beta0/beta)*dp1/pz, i.e. the full
       time of flight and not just its second-order part.  Subtracting ds here
       (as a drift inside the element would) is wrong: it would cancel the
       first-order term that Xsuite keeps. */
    if (ds != 0.0) {
        double dp1 = 1.0 + r6[delta_];
        double pz = pxyz(dp1, r6[px_], r6[py_]);
        double NormL = ds / pz;
        r6[x_] += r6[px_] * NormL;
        r6[y_] += r6[py_] * NormL;
        r6[ct_] += NormL * dp1;
    }
}

static void mis_yrot(double *r6, double phi)
{
    /* Rotation about the y axis. Forest 10.26 -- same convention as the Yrot()
       used for the pole faces. */
    if (phi != 0.0) {
        double dp1 = 1.0 + r6[delta_];
        double x = r6[x_];
        double px = r6[px_];
        double py = r6[py_];
        double c = cos(phi);
        double s = sin(phi);
        double pz = pxyz(dp1, px, py);
        double p = c*pz - s*px;
        r6[x_] = x*pz/p;
        r6[px_] = s*pz + c*px;
        r6[y_] += x*py*s/p;
        r6[ct_] += dp1*x*s/p;
    }
}

static void mis_xrot(double *r6, double phi)
{
    /* Rotation about the x axis: the x <-> y image of mis_yrot(). */
    if (phi != 0.0) {
        double dp1 = 1.0 + r6[delta_];
        double y = r6[y_];
        double px = r6[px_];
        double py = r6[py_];
        double c = cos(phi);
        double s = sin(phi);
        double pz = pxyz(dp1, px, py);
        double p = c*pz - s*py;
        r6[y_] = y*pz/p;
        r6[py_] = s*pz + c*py;
        r6[x_] += y*px*s/p;
        r6[ct_] += dp1*y*s/p;
    }
}

static void mis_srot(double *r6, double psi)
{
    /* Rotation about the s axis (roll). Purely transverse: no path length. */
    if (psi != 0.0) {
        double c = cos(psi);
        double s = sin(psi);
        double x = r6[x_];
        double y = r6[y_];
        double px = r6[px_];
        double py = r6[py_];
        r6[x_]  =  c*x + s*y;
        r6[y_]  = -s*x + c*y;
        r6[px_] =  c*px + s*py;
        r6[py_] = -s*px + c*py;
    }
}

/* ------------------------------------------------------------------ */
/* 4x4 rigid affine helpers                                            */
/* ------------------------------------------------------------------ */

static void mis_mat_mul(const double a[4][4], const double b[4][4], double out[4][4])
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            out[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j]
                      + a[i][2]*b[2][j] + a[i][3]*b[3][j];
        }
    }
}

static void mis_mat_rigid_inv(const double m[4][4], double out[4][4])
{
    /* Inverse of a rigid affine transform: transpose the rotation block and
       apply it to minus the translation. */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) out[i][j] = m[j][i];
    }
    for (int i = 0; i < 3; i++) {
        out[i][3] = -(m[0][i]*m[0][3] + m[1][i]*m[1][3] + m[2][i]*m[2][3]);
    }
    out[3][0] = 0.0; out[3][1] = 0.0; out[3][2] = 0.0; out[3][3] = 1.0;
}

static void mis_matrix(double dx, double dy, double ds,
                       double theta, double phi, double psi, double out[4][4])
{
    double s_phi = sin(phi), c_phi = cos(phi);
    double s_theta = sin(theta), c_theta = cos(theta);
    double s_psi = sin(psi), c_psi = cos(psi);

    out[0][0] = -s_phi*s_psi*s_theta + c_psi*c_theta;
    out[0][1] = -c_psi*s_phi*s_theta - c_theta*s_psi;
    out[0][2] = c_phi*s_theta;
    out[0][3] = dx;

    out[1][0] = c_phi*s_psi;
    out[1][1] = c_phi*c_psi;
    out[1][2] = s_phi;
    out[1][3] = dy;

    out[2][0] = -c_theta*s_phi*s_psi - c_psi*s_theta;
    out[2][1] = -c_psi*c_theta*s_phi + s_psi*s_theta;
    out[2][2] = c_phi*c_theta;
    out[2][3] = ds;

    out[3][0] = 0.0; out[3][1] = 0.0; out[3][2] = 0.0; out[3][3] = 1.0;
}

static void mis_frame(double part_angle, double h, double tilt_frame, double out[4][4])
{
    /* Transform from the anchor point to the element entrance (or exit) along
       the curved design orbit. */
    double c_a = cos(part_angle), s_a = sin(part_angle);
    double c_t = cos(tilt_frame), s_t = sin(tilt_frame);

    out[0][0] = (c_a - 1.0)*c_t*c_t + 1.0;
    out[0][1] = (c_a - 1.0)*c_t*s_t;
    out[0][2] = -c_t*s_a;
    out[0][3] = (c_a - 1.0)*c_t/h;

    out[1][0] = (c_a - 1.0)*c_t*s_t;
    out[1][1] = (c_a - 1.0)*s_t*s_t + 1.0;
    out[1][2] = -s_a*s_t;
    out[1][3] = (c_a - 1.0)*s_t/h;

    out[2][0] = c_t*s_a;
    out[2][1] = s_a*s_t;
    out[2][2] = c_a;
    out[2][3] = s_a/h;

    out[3][0] = 0.0; out[3][1] = 0.0; out[3][2] = 0.0; out[3][3] = 1.0;
}

/* ------------------------------------------------------------------ */
/* Entry / exit transformations                                        */
/* ------------------------------------------------------------------ */

static void exact_misalign_entry(double *r6, double dx, double dy, double ds,
        double theta, double phi, double psi, double anchor, double length,
        double angle, double h, double psi_frame)
{
    if (angle == 0.0 && (length != 0.0 || h == 0.0)) {
        /* Straight element: the anchor correction is a pure translation. */
        double mis_x = dx - anchor*cos(phi)*sin(theta);
        double mis_y = dy - anchor*sin(phi);
        double mis_s = ds - anchor*(cos(phi)*cos(theta) - 1.0);

        mis_xy_shift(r6, mis_x, mis_y);
        mis_s_shift(r6, mis_s);
        mis_yrot(r6, -theta);
        mis_xrot(r6, -phi);
        mis_srot(r6, psi);
        mis_srot(r6, psi_frame);
        return;
    }

    /* Curved element:
           misaligned_entry = F * M * inv(F)
       with M the misalignment and F the arc from the entrance to the anchor. */
    if (length != 0.0) h = angle/length;

    double mmat[4][4], fmat[4][4], finv[4][4], tmp[4][4], entry[4][4];
    mis_matrix(dx, dy, ds, theta, phi, psi, mmat);
    mis_frame(anchor*h, h, psi_frame, fmat);
    mis_mat_rigid_inv(fmat, finv);
    mis_mat_mul(fmat, mmat, tmp);
    mis_mat_mul(tmp, finv, entry);

    double mis_x = entry[0][3];
    double mis_y = entry[1][3];
    double mis_s = entry[2][3];
    double rot_theta = atan2(entry[0][2], entry[2][2]);
    double rot_phi = atan2(entry[1][2],
                           sqrt(entry[1][0]*entry[1][0] + entry[1][1]*entry[1][1]));
    double rot_psi = atan2(entry[1][0], entry[1][1]);

    mis_xy_shift(r6, mis_x, mis_y);
    mis_s_shift(r6, mis_s);
    mis_yrot(r6, -rot_theta);
    mis_xrot(r6, -rot_phi);
    mis_srot(r6, rot_psi);
    mis_srot(r6, psi_frame);
}

static void exact_misalign_exit(double *r6, double dx, double dy, double ds,
        double theta, double phi, double psi, double anchor, double length,
        double angle, double h, double psi_frame)
{
    if (angle == 0.0 && (length != 0.0 || h == 0.0)) {
        double neg = anchor - length;
        double mis_x = neg*cos(phi)*sin(theta) - dx;
        double mis_y = neg*sin(phi) - dy;
        double mis_s = neg*(cos(phi)*cos(theta) - 1.0) - ds;

        mis_srot(r6, -psi_frame);
        mis_srot(r6, -psi);
        mis_xrot(r6, phi);
        mis_yrot(r6, theta);
        mis_s_shift(r6, mis_s);
        mis_xy_shift(r6, mis_x, mis_y);
        return;
    }

    /* Curved element:
           realign = inv(F) * inv(M) * F
       with F the arc from the anchor to the element exit. */
    if (length != 0.0) h = angle/length;

    double mmat[4][4], minv[4][4], fmat[4][4], finv[4][4], tmp[4][4], realign[4][4];
    mis_matrix(dx, dy, ds, theta, phi, psi, mmat);
    mis_mat_rigid_inv(mmat, minv);
    mis_frame(angle - h*anchor, h, psi_frame, fmat);
    mis_mat_rigid_inv(fmat, finv);
    mis_mat_mul(finv, minv, tmp);
    mis_mat_mul(tmp, fmat, realign);

    double mis_x = realign[0][3];
    double mis_y = realign[1][3];
    double mis_s = realign[2][3];
    double rot_theta = atan2(realign[0][2], realign[2][2]);
    double rot_phi = atan2(realign[1][2],
                           sqrt(realign[1][0]*realign[1][0] + realign[1][1]*realign[1][1]));
    double rot_psi = atan2(realign[1][0], realign[1][1]);

    mis_srot(r6, -psi_frame);
    mis_xy_shift(r6, mis_x, mis_y);
    mis_s_shift(r6, mis_s);
    mis_yrot(r6, -rot_theta);
    mis_xrot(r6, -rot_phi);
    mis_srot(r6, rot_psi);
}

/* ------------------------------------------------------------------ */
/* Element-level helper                                                */
/* ------------------------------------------------------------------ */

struct exact_misalign {
    int active;
    double dx, dy, dz;
    double theta;   /* yaw   */
    double phi;     /* -pitch */
    double psi;     /* tilt  */
    double psi_frame;
    double anchor;
};

/* Read the standard AT geometry attributes.  Returns a filled-in descriptor
   whose "active" flag says whether an exact misalignment has to be applied. */
#define GET_EXACT_MISALIGN(ElemData, mis, Length, ExactMisalign)                 \
    do {                                                                        \
        (mis).active = 0;                                                       \
        if (ExactMisalign) {                                                    \
            (mis).dx = atGetOptionalDouble(ElemData, "dx", 0.0);                \
            (mis).dy = atGetOptionalDouble(ElemData, "dy", 0.0);                \
            (mis).dz = atGetOptionalDouble(ElemData, "dz", 0.0);                \
            (mis).theta = atGetOptionalDouble(ElemData, "yaw", 0.0);            \
            (mis).phi = -atGetOptionalDouble(ElemData, "pitch", 0.0);           \
            (mis).psi = atGetOptionalDouble(ElemData, "tilt", 0.0);             \
            (mis).psi_frame = atGetOptionalDouble(ElemData, "tilt_frame", 0.0); \
            (mis).anchor = atGetOptionalDouble(ElemData, "MisalignAnchor",      \
                                               0.5*(Length));                   \
            (mis).active = ((mis).dx != 0.0 || (mis).dy != 0.0                  \
                         || (mis).dz != 0.0 || (mis).theta != 0.0               \
                         || (mis).phi != 0.0 || (mis).psi != 0.0                \
                         || (mis).psi_frame != 0.0);                            \
        }                                                                       \
    } while (0)

#endif /*EXACT_MISALIGN_H*/
