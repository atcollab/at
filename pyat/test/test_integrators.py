"""Simple tests to make sure integrators load in Python."""
import numpy
import pytest

# noinspection PyUnresolvedReferences,PyProtectedMember
from at.lattice import Element, elements
from at import shift_elem, tilt_elem
from at import element_track, lattice_track
from at import lattice_pass, internal_lpass
from at import element_pass, internal_epass


@pytest.mark.parametrize("func", (element_track, element_pass, internal_epass))
@pytest.mark.parametrize("passmethod", ("GWigSymplecticPass", "GWigSymplecticRadPass"))
def test_gwig_symplectic_pass(rin, passmethod, func):
    # Parameters copied from one of the Diamond wigglers.
    wiggler = elements.Wiggler("w", 1.15, 0.05, 0.8)
    wiggler.PassMethod = passmethod
    func(wiggler, rin, energy=1e9)


@pytest.mark.parametrize("func", (element_track, element_pass, internal_epass))
def test_bndstrmpole_symplectic_4_pass(rin, func):
    bend = elements.Dipole("b", 1.0)
    bend.PassMethod = "BndStrMPoleSymplectic4Pass"
    func(bend, rin)


@pytest.mark.parametrize("func", (element_track, element_pass, internal_epass))
def test_pydrift(func):
    pydrift = elements.Drift("drift", 1.0, PassMethod="pyDriftPass")
    cdrift = elements.Drift("drift", 1.0, PassMethod="DriftPass")
    pyout, *_ = func(pydrift, numpy.zeros(6) + 1.0e-6)
    cout, *_ = func(cdrift, numpy.zeros(6) + 1.0e-6)
    numpy.testing.assert_equal(pyout, cout)

    shift_elem(pydrift, 1.0e-3, 1.0e-3)
    shift_elem(cdrift, 1.0e-3, 1.0e-3)
    pyout, *_ = func(pydrift, numpy.zeros(6) + 1.0e-6)
    cout, *_ = func(cdrift, numpy.zeros(6) + 1.0e-6)
    numpy.testing.assert_equal(pyout, cout)

    tilt_elem(pydrift, 1.0e-3, 1.0e-3)
    tilt_elem(cdrift, 1.0e-3, 1.0e-3)
    pyout, *_ = func(pydrift, numpy.zeros(6) + 1.0e-6)
    cout, *_ = func(cdrift, numpy.zeros(6) + 1.0e-6)
    numpy.testing.assert_equal(pyout, cout)

    # Multiple particles
    pyout, *_ = func(pydrift, numpy.zeros((6, 2)) + 1.0e-6)
    cout, *_ = func(cdrift, numpy.zeros((6, 2)) + 1.0e-6)
    numpy.testing.assert_equal(pyout, cout)


@pytest.mark.parametrize("func", (lattice_track, lattice_pass, internal_lpass))
def test_pyintegrator(hmba_lattice, func):
    params = {
        "Length": 0,
        "PassMethod": "pyIdentityPass",
    }
    id_elem = Element("py_id", **params)
    pin = numpy.zeros((6, 2)) + 1.0e-6
    pout1, *_ = func(hmba_lattice, pin.copy(), nturns=1)
    hmba_lattice = hmba_lattice + [id_elem]
    pout2, *_ = func(hmba_lattice, pin.copy(), nturns=1)
    numpy.testing.assert_equal(pout1, pout2)
