import math

import numpy as np
from numpy.testing import assert_allclose

import at
from at.lattice import get_geometry3


def _build_3d_lattice():
    """Build a 3D lattice."""

    def set_R1_R2(elem: at.Element, rots: float) -> None:
        """Set the R1 and R2 roation matrices."""
        cs = np.cos(rots)
        sn = np.sin(rots)
        rm = np.asfortranarray(np.diag([cs, cs, cs, cs, 1.0, 1.0]))
        rm[0, 2] = sn
        rm[1, 3] = sn
        rm[2, 0] = -sn
        rm[3, 1] = -sn
        elem.R1 = rm
        elem.R2 = rm.T

    dip = at.Dipole("dip", 1.0, math.pi / 4.0)
    dr1 = at.Drift("d1", 2.0)
    dr2 = at.Drift("d2", 1.0)
    dup = at.Dipole("dup", 1.0, math.pi / 32.0)
    set_R1_R2(dup, -math.pi / 2.0)
    ddn = at.Dipole("ddn", 1.0, math.pi / 32.0)
    set_R1_R2(ddn, math.pi / 2.0)
    vup = [dup, dr2, ddn, dip]
    vdn = [ddn, dr2, dup, dip]
    return at.Lattice(
        [dr2, dip] + [dr1, dip] + vup + [dr1, dip] * 3 + vdn + [dr1, dip] + [dr2],
        energy=1.0e9,
    )


def test_geometry():
    expected = np.rec.array(
        [
            (0.00000000e00, 0.00000000e00, 0.00000000e00, 0.00000000e00),
            (1.00000000e00, 0.00000000e00, 0.00000000e00, 0.00000000e00),
            (1.90031632e00, -3.72923229e-01, 0.00000000e00, -7.85398163e-01),
            (3.31452988e00, -1.78713679e00, 0.00000000e00, -7.85398163e-01),
            (3.68745311e00, -2.68745311e00, 0.00000000e00, -1.57079633e00),
            (3.68745311e00, -3.68584750e00, 4.90479714e-02, -1.57079633e00),
            (3.68745311e00, -4.68103223e00, 1.47065112e-01, -1.57079633e00),
            (3.68745311e00, -5.67942662e00, 1.96113083e-01, -1.57079633e00),
            (3.31452988e00, -6.57974294e00, 1.96113083e-01, -2.35619449e00),
            (1.90031632e00, -7.99395650e00, 1.96113083e-01, -2.35619449e00),
            (1.00000000e00, -8.36687973e00, 1.96113083e-01, 3.14159265e00),
            (-1.00000000e00, -8.36687973e00, 1.96113083e-01, 3.14159265e00),
            (-1.90031632e00, -7.99395650e00, 1.96113083e-01, 2.35619449e00),
            (-3.31452988e00, -6.57974294e00, 1.96113083e-01, 2.35619449e00),
            (-3.68745311e00, -5.67942662e00, 1.96113083e-01, 1.57079633e00),
            (-3.68745311e00, -4.68103223e00, 1.47065112e-01, 1.57079633e00),
            (-3.68745311e00, -3.68584750e00, 4.90479714e-02, 1.57079633e00),
            (-3.68745311e00, -2.68745311e00, 0.00000000e00, 1.57079633e00),
            (-3.31452988e00, -1.78713679e00, 7.10533554e-18, 7.85398163e-01),
            (-1.90031632e00, -3.72923229e-01, 1.92434190e-17, 7.85398163e-01),
            (-1.00000000e00, -7.77156117e-16, 2.30661980e-17, -3.66943061e-16),
            (2.99760217e-15, -1.14409918e-15, 2.42455998e-17, -3.66943061e-16),
        ],
        dtype=[("x", "<f8"), ("y", "<f8"), ("z", "<f8"), ("angle", "<f8")],
    )

    ring = _build_3d_lattice()
    g2 = get_geometry3(ring)
    assert_allclose(g2.x, expected.x, atol=1.0e-12)
    assert_allclose(g2.y, expected.y, atol=1.0e-12)
    assert_allclose(g2.z, expected.z, atol=1.0e-12)
    assert_allclose(g2.angle, expected.angle, atol=1.0e-12)
