"""Exact misalignment of the Exact* pass methods, checked against Xsuite.

These tests are skipped unless xtrack is installed.  They convert the AT
element with at.line_from_lattice and require the two codes to agree to a
tolerance well below what AT's linearised R1/T1 misalignment achieves.
"""

import numpy as np
import pytest

import at

xt = pytest.importorskip("xtrack")

ENERGY = 10e9
TOL = 1e-9  # AT's linearised misalignment is 1e-4..1e-3 on these cases


def _particles(amp=5e-3, n=7):
    a = np.linspace(-amp, amp, n)
    g1, g2 = np.meshgrid(a, a)
    p = np.zeros((6, g1.size))
    p[0] = g1.ravel()
    p[2] = g2.ravel()
    return p


def _elements():
    yield "multipole", at.Multipole(
        "M", 1.0, [0.0, 0.0, 0.0, 0.0], [0.0, 0.5, 0.0, 0.0],
        PassMethod="ExactMultipolePass", NumIntSteps=200)
    yield "sector", at.Dipole(
        "D", 2.0, 1.0, k=1.0, PassMethod="ExactSectorBendPass",
        NumIntSteps=200, EntranceAngle=0.5, ExitAngle=0.5,
        FringeBendEntrance=1, FringeBendExit=1,
        FringeQuadEntrance=1, FringeQuadExit=1)
    yield "rectangular", at.Dipole(
        "R", 2.0, 0.4, k=0.3, PassMethod="ExactRectangularBendPass",
        NumIntSteps=200, EntranceAngle=0.2, ExitAngle=0.2,
        FringeBendEntrance=1, FringeBendExit=1,
        FringeQuadEntrance=1, FringeQuadExit=1)


MISALIGNMENTS = [
    dict(dx=1e-3),
    dict(dy=1e-3),
    dict(dz=1e-3),
    dict(tilt=1e-3),
    dict(pitch=1e-3),
    dict(yaw=1e-3),
    dict(dx=2e-3, dy=-1e-3, dz=1e-3, tilt=1e-2, pitch=-5e-3, yaw=8e-3),
]


def _track_both(elem):
    lattice = at.Lattice([elem], energy=ENERGY, periodicity=1)
    line = at.line_from_lattice(lattice, match_model=True)
    if line.particle_ref is None:
        line.particle_ref = xt.Particles(p0c=ENERGY, mass0=xt.ELECTRON_MASS_EV)

    parts = _particles()
    pout, *_ = lattice.track(np.asfortranarray(parts.copy()))
    at_out = np.asarray(pout).reshape(6, parts.shape[1], -1)[:, :, -1]

    p = line.build_particles(
        x=parts[0], px=parts[1], y=parts[2], py=parts[3],
        delta=parts[4], zeta=parts[5])
    line.track(p)
    keep = p.state > 0
    idx = p.particle_id[keep]
    xs_out = np.full_like(at_out, np.nan)
    for i, attr in enumerate(["x", "px", "y", "py", "delta"]):
        xs_out[i, idx] = getattr(p, attr)[keep]
    xs_out[5, idx] = -p.zeta[keep]
    return at_out, xs_out


@pytest.mark.parametrize("misalign", MISALIGNMENTS)
@pytest.mark.parametrize("name,elem", list(_elements()), ids=lambda e: getattr(e, "FamName", e))
def test_exact_misalignment_matches_xsuite(name, elem, misalign):
    elem = elem.copy()
    at.transform_elem(elem, **misalign)
    elem.ExactMisalign = 1
    at_out, xs_out = _track_both(elem)
    finite = np.isfinite(at_out - xs_out).all(axis=0)
    assert finite.any(), "all particles lost"
    assert np.nanmax(np.abs((at_out - xs_out)[:, finite])) < TOL


def test_k0h_correction_matches_xsuite():
    """A dipole field error in a curved frame needs its (1 + h x) correction."""
    elem = at.Dipole("D", 2.0, 1.0, k=1.0, PassMethod="ExactSectorBendPass",
                     NumIntSteps=200, EntranceAngle=0.5, ExitAngle=0.5,
                     FringeBendEntrance=1, FringeBendExit=1,
                     FringeQuadEntrance=1, FringeQuadExit=1)
    elem.PolynomB[0] += 1e-4 * abs(elem.K)
    at_out, xs_out = _track_both(elem)
    finite = np.isfinite(at_out - xs_out).all(axis=0)
    assert np.nanmax(np.abs((at_out - xs_out)[:, finite])) < TOL
