"""Simple tests to make sure integrators load in Python."""
import numpy
import pytest

from at import atpass, elements
from at.lattice import Lattice
from at import Element, lattice_pass
from at import set_shift, set_tilt


def test_exact_hamiltonian_pass(rin):
    drift = elements.Multipole('m1', 1, [0, 0, 0, 0], [0, 0, 0, 0])
    drift.Type = 0
    drift.PassMethod = 'ExactHamiltonianPass'
    drift.BendingAngle = 0
    l = Lattice([drift], name='lat', energy=3e9)
    atpass(l, rin, 1)


def test_exact_hamiltonian_pass_with_dls_dipole(rin):
    bend = elements.Multipole('rb', 0.15, [0, 0, 0, 0], [-0.0116333, 3.786786, 0, 0])
    bend.Type = 1
    bend.PassMethod = 'ExactHamiltonianPass'
    bend.BendingAngle = -0.001745
    bend.Energy = 3.5e9
    bend.MaxOrder = 3
    l = Lattice([bend], name='lat', energy=3.5e9)
    atpass(l, rin, 1)
    # Results from Matlab
    expected = numpy.array([9.23965e-9, 1.22319e-5, 0, 0, 0, -4.8100e-10]).reshape(6,1)
    numpy.testing.assert_allclose(rin, expected, rtol=1e-5, atol=1e-6)


@pytest.mark.parametrize('passmethod', ('GWigSymplecticPass', 'GWigSymplecticRadPass'))
def test_gwig_symplectic_pass(rin, passmethod):
    # Parameters copied from one of the Diamond wigglers.
    wiggler = elements.Wiggler('w', 1.15, 0.05, 0.8, 3e9)
    wiggler.PassMethod = passmethod
    l = Lattice([wiggler], name='lat', energy=3e9)
    atpass(l, rin, 1)


def test_bndstrmpole_symplectic_4_pass(rin):
    bend = elements.Dipole('b', 1.0)
    bend.PassMethod = 'BndStrMPoleSymplectic4Pass'
    l = Lattice([bend], name='lat', energy=3e9)
    atpass(l, rin, 1)


def test_pydrift(rin):
    pydrift = [elements.Drift('drift', 1.0, 
               PassMethod='pyDriftPass')]
    cdrift = [elements.Drift('drift', 1.0,
               PassMethod='DriftPass')]
    pyout = lattice_pass(pydrift, rin, nturns=1)
    cout = lattice_pass(cdrift, rin, nturns=1)
    numpy.testing.assert_equal(pyout, cout)

    set_shift(pydrift, [1.0e-3], [1.0e-3], relative=False)
    set_shift(cdrift, [1.0e-3], [1.0e-3], relative=False)
    pyout = lattice_pass(pydrift, rin, nturns=1)
    cout = lattice_pass(cdrift, rin, nturns=1)
    numpy.testing.assert_equal(pyout, cout)

    set_tilt(pydrift, [1.0e-3], relative=False)
    set_tilt(cdrift, [1.0e-3], relative=False)
    pyout = lattice_pass(pydrift, rin, nturns=1)
    cout = lattice_pass(cdrift, rin, nturns=1)
    numpy.testing.assert_equal(pyout, cout)


def test_pyintegrator(hmba_lattice):
    params = {'Length': 0,
              'PassMethod': 'pyIdentityPass',
              }
    id_elem = Element('py_id', **params)
    pin = numpy.zeros(6)+1.0e-6
    pout1 = lattice_pass(hmba_lattice, pin.copy(), nturns=1)
    pin = numpy.zeros(6)+1.0e-6
    hmba_lattice = hmba_lattice + [id_elem]
    pout2 = lattice_pass(hmba_lattice, pin, nturns=1)
    numpy.testing.assert_equal(pout1, pout2)
