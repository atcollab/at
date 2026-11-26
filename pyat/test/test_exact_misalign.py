"""Exact misalignment of the magnet pass methods (ExactMisalign=1).

Reference: M. Venturini and T. Hellert, "Particle-ray tracing through misaligned
lattice elements in accelerators", https://doi.org/10.1016/j.nima.2022.167487
"""

import numpy as np
import pytest

import at
from at import ReferencePoint
from at.tracking import diffusion_matrix

MISALIGNMENTS = {
    "dx": {"dx": 1e-3},
    "dy": {"dy": -1e-3},
    "dz": {"dz": 1e-3},
    "tilt": {"tilt": 1e-3},
    "pitch": {"pitch": 1e-3},
    "yaw": {"yaw": -1e-3},
    "all": {"dx": 2e-3, "dy": -1e-3, "dz": 1e-3, "tilt": 10e-3, "pitch": -5e-3,
            "yaw": 8e-3},
}


def sector_bend(passmethod="ExactSectorBendPass"):
    return at.Dipole("B", 2.0, 0.5, 0.3, PassMethod=passmethod, NumIntSteps=50,
                     PolynomB=[0.0, 0.3, 1.0], EntranceAngle=0.2, ExitAngle=0.1,
                     Energy=6e9)


def rect_bend(passmethod="ExactRectangularBendPass"):
    b = at.Dipole("R", 2.0, 0.4, PassMethod=passmethod,
                  NumIntSteps=50, EntranceAngle=0.2, ExitAngle=0.2, Energy=6e9)
    b.rbendtune()
    return b


def multipole(passmethod="ExactMultipolePass"):
    return at.Multipole("M", 1.0, [0.0, 1e-3, 0.0], [0.0, 0.5, 2.0],
                        PassMethod=passmethod, NumIntSteps=50, Energy=6e9)


def _factory(make, passmethod):
    return lambda: make(passmethod)


ELEMENTS = {
    pm: _factory(make, pm)
    for make, pms in (
        (multipole, ["ExactMultipolePass", "ExactMultipoleRadPass",
                     "StrMPoleSymplectic4Pass", "StrMPoleSymplectic4RadPass",
                     "StrMPoleSymplectic6Pass", "StrMPoleSymplectic6RadPass"]),
        (sector_bend, ["ExactSectorBendPass", "ExactSectorBendRadPass",
                       "BndMPoleSymplectic4Pass", "BndMPoleSymplectic4RadPass",
                       "BndMPoleSymplectic6Pass", "BndMPoleSymplectic6RadPass",
                       "BndMPoleSymplectic4E2Pass", "BndMPoleSymplectic4E2RadPass"]),
        (rect_bend, ["ExactRectangularBendPass", "ExactRectangularBendRadPass",
                     "ExactRectBendPass", "ExactRectBendRadPass"]),
    )
    for pm in pms
}
# Pass methods propagating the diffusion matrix
DIFFUSION = ["ExactMultipoleRadPass", "StrMPoleSymplectic4RadPass",
             "StrMPoleSymplectic6RadPass", "BndMPoleSymplectic4RadPass",
             "BndMPoleSymplectic6RadPass"]


def track(elem, r):
    """Track one 6-vector or an array of 6-vectors through one element."""
    r = np.asfortranarray(np.atleast_2d(np.asarray(r, dtype=float)).T.copy())
    return elem.track(r, energy=getattr(elem, "Energy", 6e9)).reshape(6, -1).T


def misaligned(factory, reference=ReferencePoint.CENTRE, exact=True, **kw):
    elem = factory()
    at.transform_elem(elem, reference=reference, **kw)
    if exact:
        elem.ExactMisalign = 1
    return elem


# --- Independent implementation of the paper (Eqs. 3-7, 12, 18-26, 31-37) ---


def _ru(ax, ay, az):
    """Eq. (3): R_u = Rx(ax) Ry(ay) Rz(az)."""
    cx, sx, cy, sy, cz, sz = np.cos(ax), np.sin(ax), np.cos(ay), np.sin(ay), \
        np.cos(az), np.sin(az)
    return np.array([
        [cy * cz, -cy * sz, sy],
        [sx * sy * cz + cx * sz, cx * cz - sx * sy * sz, -sx * cy],
        [sx * sz - cx * sy * cz, cx * sy * sz + sx * cz, cx * cy],
    ])


def _rb(theta):
    """Eq. (12) without tilt."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])


def paper_geometry(length, angle, reference, dx=0.0, dy=0.0, dz=0.0, tilt=0.0,
                   pitch=0.0, yaw=0.0):
    """Rotation and displacement of the entrance and exit faces."""
    ra = _ru(pitch, yaw, tilt)
    da = np.array([dx, dy, dz])
    z = np.array([0.0, 0.0, 1.0])
    if reference is ReferencePoint.ENTRANCE:  # Sec. 2
        r, d = ra, da
    else:  # Sec. 5: origin at the midpoint of the chord
        rh = _rb(angle / 2)
        r = rh @ ra @ rh.T  # Eq. (31)
        half_chord = length / 2 if angle == 0 else length / angle * np.sin(angle / 2)
        d = half_chord * rh @ z + rh @ da - half_chord * r @ rh @ z  # Eqs. (33-36)
    rb = _rb(angle)
    if angle == 0:
        opp = np.array([0.0, 0.0, length])  # Eq. (26)
    else:
        rc = length / angle
        opp = np.array([rc * (np.cos(angle) - 1), 0.0, rc * np.sin(angle)])  # Eq. (24)
    r_exit = rb.T @ r.T @ rb  # Eq. (18)
    d_exit = (r @ rb).T @ (opp - (r @ opp + d))  # Eqs. (23, 25)
    return (r, d), (r_exit, d_exit)


def paper_map(r6, rot, d):
    """Eqs. (4-7), ultra-relativistic limit (Appendix A)."""
    x, px, y, py, delta, ct = r6
    pz = np.sqrt((1 + delta) ** 2 - px**2 - py**2)
    rr = rot.T @ (np.array([x, y, 0.0]) - d)
    pp = rot.T @ np.array([px, py, pz])
    lg = -rr[2]
    return np.array([rr[0] + lg * pp[0] / pp[2], pp[0], rr[1] + lg * pp[1] / pp[2],
                     pp[1], delta, ct + lg * (1 + delta) / pp[2]])


def _geometry(elem, reference, mis):
    return paper_geometry(elem.Length, getattr(elem, "BendingAngle", 0.0),
                          reference, **mis)


PARTICLES = np.array([
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1e-3, 1e-3, -1e-3, 5e-4, 0.0, 0.0],
    [-2e-3, 5e-4, 1e-3, -1e-3, 1e-2, 1e-3],
])


# --- Tests -------------------------------------------------------------------


@pytest.mark.parametrize("name", ELEMENTS)
def test_no_misalignment_identical(name):
    """With zero misalignment, the exact maps are the identity (to round-off
    for bends, where R = R_B/2 R_B/2^T)."""
    ref = track(ELEMENTS[name](), PARTICLES)
    out = track(misaligned(ELEMENTS[name]), PARTICLES)
    np.testing.assert_allclose(out, ref, rtol=0, atol=2e-15)


@pytest.mark.parametrize("name", ELEMENTS)
def test_flag_zero_is_linear(name):
    """ExactMisalign=0 is identical to no flag."""
    elem = misaligned(ELEMENTS[name], **MISALIGNMENTS["all"])
    elem.ExactMisalign = 0
    lin = misaligned(ELEMENTS[name], exact=False, **MISALIGNMENTS["all"])
    np.testing.assert_array_equal(track(elem, PARTICLES), track(lin, PARTICLES))


@pytest.mark.parametrize("reference", [ReferencePoint.CENTRE, ReferencePoint.ENTRANCE])
@pytest.mark.parametrize("mis", MISALIGNMENTS)
@pytest.mark.parametrize("name", ELEMENTS)
def test_matches_paper(name, mis, reference):
    """Exact tracking == paper exit map o aligned element o paper entrance map."""
    elem = misaligned(ELEMENTS[name], reference, **MISALIGNMENTS[mis])
    aligned = ELEMENTS[name]()
    (r1, d1), (r2, d2) = _geometry(aligned, reference, MISALIGNMENTS[mis])
    for r6 in PARTICLES:
        expected = paper_map(track(aligned, paper_map(r6, r1, d1))[0], r2, d2)
        np.testing.assert_allclose(track(elem, r6)[0], expected, rtol=0, atol=1e-14)


@pytest.mark.parametrize("tilt_frame", [0.0, 0.3])
@pytest.mark.parametrize("name", ["ExactMultipolePass", "ExactSectorBendPass"])
def test_tilt_frame(name, tilt_frame):
    """A frame tilt alone is a rotation about Z: exact == linear."""
    elem = misaligned(ELEMENTS[name], tilt_frame=tilt_frame)
    lin = misaligned(ELEMENTS[name], exact=False, tilt_frame=tilt_frame)
    np.testing.assert_allclose(track(elem, PARTICLES), track(lin, PARTICLES),
                               rtol=0, atol=1e-15)


@pytest.mark.parametrize("mis", MISALIGNMENTS)
def test_misaligned_drift_is_a_drift(mis):
    """A field-free element is a drift whatever its misalignment (paper Sec. 6).

    Exact with the exact maps, approximate with the linear ones.
    """
    def drift():
        return at.Multipole("D", 1.5, [0.0], [0.0], PassMethod="ExactMultipolePass",
                            NumIntSteps=1)

    ref = track(drift(), PARTICLES)
    exact = track(misaligned(drift, **MISALIGNMENTS[mis]), PARTICLES)
    linear = track(misaligned(drift, exact=False, **MISALIGNMENTS[mis]), PARTICLES)
    np.testing.assert_allclose(exact, ref, rtol=0, atol=1e-15)
    if mis in ("tilt", "dx", "dy", "dz"):  # linear maps are exact for these
        np.testing.assert_allclose(linear, ref, rtol=0, atol=1e-15)
    else:
        assert np.abs(linear - ref).max() > 1e-10


@pytest.mark.parametrize("tilt_frame", [0.0, 0.3])
@pytest.mark.parametrize("reference", [ReferencePoint.CENTRE, ReferencePoint.ENTRANCE])
@pytest.mark.parametrize("name", ["ExactMultipolePass", "ExactSectorBendPass",
                                  "ExactRectangularBendPass", "StrMPoleSymplectic6Pass",
                                  "BndMPoleSymplectic6Pass"])
def test_linear_limit(name, reference, tilt_frame):
    """Exact and linear maps share the geometry: their difference is of higher order.

    The linear maps are wrong at second order in the coordinates at each face,
    with coefficients proportional to the misalignment. Inside the element the
    coordinates are of the order of the misalignment itself, so scaling both the
    misalignment and the particle amplitude by k scales the difference by k^3.
    """
    direction = np.array([1.0, 0.5, -1.0, 0.7, 0.0, 0.0])
    diffs = []
    for k in (1.0, 0.1):
        kw = {key: k * v for key, v in MISALIGNMENTS["all"].items()}
        exact = misaligned(ELEMENTS[name], reference, tilt_frame=tilt_frame, **kw)
        linear = misaligned(ELEMENTS[name], reference, exact=False,
                            tilt_frame=tilt_frame, **kw)
        r6 = k * 1e-3 * direction
        diffs.append(np.abs(track(exact, r6) - track(linear, r6)).max())
    assert 700 < diffs[0] / diffs[1] < 1300


@pytest.mark.parametrize("name", ["ExactMultipolePass", "ExactSectorBendPass",
                                  "ExactRectangularBendPass", "StrMPoleSymplectic6Pass",
                                  "BndMPoleSymplectic6Pass"])
def test_symplectic(name):
    elem = misaligned(ELEMENTS[name], **MISALIGNMENTS["all"])
    r0 = np.array([1e-3, -5e-4, 5e-4, 1e-3, 1e-3, 0.0])
    eps = 1e-6
    m = np.zeros((6, 6))
    for j in range(6):
        dr = np.zeros(6)
        dr[j] = eps
        m[:, j] = (track(elem, r0 + dr)[0] - track(elem, r0 - dr)[0]) / (2 * eps)
    s = at.jmat(3)
    np.testing.assert_allclose(m.T @ s @ m, s, rtol=0, atol=1e-8)


@pytest.mark.parametrize("name", DIFFUSION)
def test_diffusion_matrix(name):
    """The diffusion matrix goes through the Jacobian of the exact exit map.

    Expected: J B_body J^T, with B_body the diffusion matrix of the aligned
    element along the orbit inside the misaligned one, and J the numerical
    Jacobian of the paper's exit map at the body exit.
    """
    mis = MISALIGNMENTS["all"]
    elem = misaligned(ELEMENTS[name], **mis)
    aligned = ELEMENTS[name]()
    (r1, d1), (r2, d2) = _geometry(aligned, ReferencePoint.CENTRE, mis)
    orbit = np.array([1e-3, -2e-4, 5e-4, 3e-4, 1e-3, 0.0])
    inside = paper_map(orbit, r1, d1)
    body_exit = track(aligned, inside)[0]
    b_body = diffusion_matrix(aligned, inside, energy=aligned.Energy)
    eps = 1e-7
    jac = np.zeros((6, 6))
    for j in range(6):
        dr = np.zeros(6)
        dr[j] = eps
        jac[:, j] = (paper_map(body_exit + dr, r2, d2)
                     - paper_map(body_exit - dr, r2, d2)) / (2 * eps)
    expected = jac @ b_body @ jac.T
    got = diffusion_matrix(elem, orbit, energy=elem.Energy)
    np.testing.assert_allclose(got, expected, rtol=1e-7, atol=1e-8 * np.abs(expected).max())
    # The linear R2 gives a different result at this misalignment
    lin = diffusion_matrix(misaligned(ELEMENTS[name], exact=False, **mis), orbit,
                           energy=elem.Energy)
    assert np.abs(lin - got).max() > 1e-9 * np.abs(got).max()


def test_non_rigid_matrix_rejected():
    """R1 that is not a rotation cannot be interpreted as a misalignment."""
    elem = multipole()
    elem.R1 = 2.0 * np.eye(6)
    elem.ExactMisalign = 1
    with pytest.raises(ValueError, match="ExactMisalign"):
        track(elem, PARTICLES)


def test_divide():
    """Slices keep the flag; inner slices have no transformation and are unchanged."""
    elem = misaligned(multipole, **MISALIGNMENTS["all"])
    parts = elem.divide([0.3, 0.3, 0.4])
    assert all(p.ExactMisalign == 1 for p in parts)
    inner = parts[1]
    ref = inner.deepcopy()
    del ref.ExactMisalign
    np.testing.assert_array_equal(track(inner, PARTICLES), track(ref, PARTICLES))
    whole = track(elem, PARTICLES)
    sliced = PARTICLES
    for p in parts:
        sliced = track(p, sliced)
    np.testing.assert_allclose(sliced, whole, rtol=0, atol=1e-9)
