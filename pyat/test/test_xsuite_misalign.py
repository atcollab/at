"""Misalignments in the Xsuite converter (no xtrack needed)."""

from __future__ import annotations

import warnings

import at
import numpy as np
import pytest
from at.load.xsuite_elements import XsElement
from numpy.testing import assert_allclose

LARGE = dict(dx=2e-3, dy=-1e-3, dz=1e-3, tilt=1e-2, pitch=-5e-3, yaw=8e-3)


def _quad():
    return at.Quadrupole("Q", 1.0, 0.5, PassMethod="ExactMultipolePass")


def _sbend():
    return at.Dipole("B", 2.0, 1.0, PassMethod="ExactSectorBendPass")


def _misaligned(make, reference, **kwargs):
    elem = make()
    at.transform_elem(elem, reference=reference, **kwargs)
    return elem


@pytest.mark.parametrize("make", [_quad, _sbend])
@pytest.mark.parametrize(
    "reference", [at.ReferencePoint.CENTRE, at.ReferencePoint.ENTRANCE]
)
@pytest.mark.parametrize(
    "mis", [dict(dx=1e-3), dict(tilt=1e-3), dict(pitch=1e-3), dict(yaw=1e-3), LARGE]
)
def test_roundtrip(make, reference, mis):
    """AT -> Xsuite -> AT gives back the same transformation matrices."""
    elem = _misaligned(make, reference, **mis)
    back = XsElement.static_from_at(elem).to_at()
    for attr in ("T1", "R1", "T2", "R2"):
        assert_allclose(getattr(back, attr), getattr(elem, attr), rtol=0, atol=1e-15)


@pytest.mark.parametrize(
    "reference", [at.ReferencePoint.CENTRE, at.ReferencePoint.ENTRANCE]
)
def test_reference_point(reference):
    """The pivot is at the entrance for ENTRANCE (exported as CENTRE before)."""
    xs = XsElement.static_from_at(_misaligned(_quad, reference, tilt=1e-3))
    expected = 0.0 if reference is at.ReferencePoint.ENTRANCE else 0.5
    assert xs["rot_shift_anchor"] == expected


@pytest.mark.parametrize(
    ("attr", "xsattr", "sign"),
    [("tilt", "rot_s_rad_no_frame", 1.0), ("pitch", "rot_x_rad", -1.0),
     ("yaw", "rot_y_rad", 1.0)],
)
def test_single_angle(attr, xsattr, sign):
    """With one angle on a straight element the two conventions coincide."""
    xs = XsElement.static_from_at(_misaligned(_quad, at.ReferencePoint.CENTRE,
                                              **{attr: 1e-3}))
    assert xs[xsattr] == pytest.approx(sign * 1e-3, abs=1e-18)
    assert all(xs.get(k, 0.0) == 0.0 for k in ("shift_x", "shift_y", "shift_s"))


def test_rotation_order():
    """Xsuite composes Ry Rx Rz, AT Rx Ry Rz: they differ at O(pitch * yaw)."""
    xs = XsElement.static_from_at(_misaligned(_quad, at.ReferencePoint.CENTRE,
                                              pitch=1e-2, yaw=1e-2))
    assert abs(xs["rot_s_rad_no_frame"]) == pytest.approx(1e-2 * 1e-2, rel=1e-3)


def test_bend_pivot():
    """A bend rotates around the chord midpoint in AT, the arc midpoint in Xsuite:
    the change of pivot is absorbed in the shifts."""
    tilt = 1e-3
    xs = XsElement.static_from_at(_misaligned(_sbend, at.ReferencePoint.CENTRE,
                                              tilt=tilt))
    sagitta = 2.0 / 1.0 * (1.0 - np.cos(0.5))
    assert xs["rot_s_rad_no_frame"] == pytest.approx(tilt, abs=1e-18)
    assert xs["shift_x"] == pytest.approx(-sagitta * (1.0 - np.cos(tilt)), rel=1e-9)
    # (I - Rz(tilt)) (-sagitta, 0, 0)
    assert xs["shift_y"] == pytest.approx(sagitta * np.sin(tilt), rel=1e-9)


@pytest.mark.parametrize(
    ("xsk0", "expected"),
    [
        (dict(k0=0.5, k0_from_h=True), 1e-3),  # as written by Xsuite: k0 = h
        (dict(k0_from_h=True), 1e-3),
        (dict(k0=0.501, k0_from_h=False), 2e-3),  # k0 is the total dipole field
        (dict(k0=0.501), 2e-3),
    ],
)
def test_bend_k0(xsk0, expected):
    """AT PolynomB[0] excludes the curvature h, Xsuite k0 is the total field."""
    xs = XsElement.static_from_dict(
        dict(__class__="Bend", length=2.0, angle=1.0, knl=[2e-3, 0, 0], **xsk0),
        name="B",
    )
    assert xs.to_at().PolynomB[0] == pytest.approx(expected, abs=1e-15)


@pytest.mark.parametrize(
    "reference", [at.ReferencePoint.CENTRE, at.ReferencePoint.ENTRANCE]
)
def test_rbend_anchor(reference):
    """Xsuite measures rot_shift_anchor along the arc, also for an RBend."""
    elem = at.Dipole("B", 1.0, 1.0, EntranceAngle=0.5, ExitAngle=0.5,
                     PassMethod="ExactRectangularBendPass")
    at.transform_elem(elem, reference=reference, tilt=1e-3)
    xs = XsElement.static_from_at(elem)
    expected = 0.0 if reference is at.ReferencePoint.ENTRANCE else 0.5
    assert type(xs).__name__ == "RBend"
    assert xs["rot_shift_anchor"] == pytest.approx(expected, abs=1e-15)
    back = xs.to_at()
    for attr in ("T1", "R1", "T2", "R2"):
        assert_allclose(getattr(back, attr), getattr(elem, attr), rtol=0, atol=1e-15)


@pytest.mark.parametrize(
    ("passmethod", "integrator", "kicks"),
    [("StrMPoleSymplectic4Pass", "yoshida-4", 30), ("StrMPoleSymplectic6Pass",
     "yoshida-6", 70), ("ExactMultipolePass", "yoshida-6", 70)],
)
def test_integrator(passmethod, integrator, kicks):
    """The Yoshida order follows the PassMethod, as NumIntKicks does."""
    elem = at.Quadrupole("Q", 1.0, 0.5, PassMethod=passmethod, NumIntSteps=10)
    xs = XsElement.static_from_at(elem)
    assert xs["integrator"] == integrator
    assert xs["num_multipole_kicks"] == kicks


@pytest.mark.parametrize(
    ("make", "passmethod"),
    [(at.Dipole, pm) for pm in (
        "BndMPoleSymplectic4Pass", "BndMPoleSymplectic6Pass", "ExactSectorBendPass",
        "ExactRectBendPass", "ExactRectangularBendPass")]
    + [(at.Quadrupole, pm) for pm in (
        "StrMPoleSymplectic4Pass", "StrMPoleSymplectic6Pass", "ExactMultipolePass")],
)
def test_passmethod_roundtrip(make, passmethod):
    """AT -> Xsuite -> AT keeps the PassMethod and NumIntSteps."""
    elem = make("E", 1.0, 0.5, PassMethod=passmethod, NumIntSteps=123)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        back = XsElement.static_from_at(elem).to_at()
    assert (back.PassMethod, back.NumIntSteps) == (passmethod, 123)


@pytest.mark.parametrize(
    ("passmethod", "model"),
    [("ExactRectangularBendPass", "drift-kick-drift-exact"),
     ("ExactRectangularBendRadPass", "drift-kick-drift-exact"),
     ("ExactRectBendPass", "bend-kick-bend"), ("ExactRectBendRadPass", "bend-kick-bend")],
)
def test_rectangular_bends(passmethod, model):
    """All rectangular bends become straight-body RBends."""
    elem = at.Dipole("B", 1.0, 1.0, EntranceAngle=0.5, ExitAngle=0.5,
                     PassMethod=passmethod)
    xs = XsElement.static_from_at(elem)
    assert type(xs).__name__ == "RBend"
    assert (xs["model"], xs["rbend_model"]) == (model, "straight-body")
