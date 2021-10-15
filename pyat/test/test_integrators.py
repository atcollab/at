"""Simple tests to make sure integrators load in Python."""
import numpy
import pytest

# noinspection PyUnresolvedReferences,PyProtectedMember
from at.tracking import _atpass as atpass
from at.tracking import lattice_pass, element_pass
from at.lattice import Lattice, Element, elements
from at import set_shift, set_tilt, shift_elem, tilt_elem


def test_exact_hamiltonian_pass(rin):
    drift = elements.Multipole('m1', 1, [0, 0, 0, 0], [0, 0, 0, 0])
    drift.Type = 0
    drift.PassMethod = 'ExactHamiltonianPass'
    drift.BendingAngle = 0
    element_pass(drift, rin)


def test_exact_hamiltonian_pass_with_dls_dipole(rin):
    bend = elements.Multipole('rb', 0.15, [0, 0, 0, 0], [-0.0116333, 3.786786, 0, 0])
    bend.Type = 1
    bend.PassMethod = 'ExactHamiltonianPass'
    bend.BendingAngle = -0.001745
    bend.Energy = 3.5e9
    bend.MaxOrder = 3
    element_pass(bend, rin)
    # Results from Matlab
    expected = numpy.array([9.23965e-9, 1.22319e-5, 0, 0, 0, -4.8100e-10]).reshape(6, 1)
    numpy.testing.assert_allclose(rin, expected, rtol=1e-5, atol=1e-6)


@pytest.mark.parametrize('passmethod', ('GWigSymplecticPass', 'GWigSymplecticRadPass'))
def test_gwig_symplectic_pass(rin, passmethod):
    # Parameters copied from one of the Diamond wigglers.
    wiggler = elements.Wiggler('w', 1.15, 0.05, 0.8, 3e9)
    wiggler.PassMethod = passmethod
    element_pass(wiggler, rin)


def test_bndstrmpole_symplectic_4_pass(rin):
    bend = elements.Dipole('b', 1.0)
    bend.PassMethod = 'BndStrMPoleSymplectic4Pass'
    element_pass(bend, rin)


def test_pydrift():
    pydrift = elements.Drift('drift', 1.0, PassMethod='pyDriftPass')
    cdrift = elements.Drift('drift', 1.0, PassMethod='DriftPass')
    pylattice = Lattice([pydrift], energy=1.0E9)
    clattice = Lattice([cdrift], energy=1.0E9)
    pyout = lattice_pass(pylattice, numpy.zeros(6)+1.0e-6, nturns=1)
    cout = lattice_pass(clattice, numpy.zeros(6)+1.0e-6, nturns=1)
    numpy.testing.assert_equal(pyout, cout)

    set_shift(pylattice, [1.0e-3], [1.0e-3], relative=False)
    set_shift(clattice, [1.0e-3], [1.0e-3], relative=False)
    pyout = lattice_pass(pylattice, numpy.zeros(6)+1.0e-6, nturns=1)
    cout = lattice_pass(clattice, numpy.zeros(6)+1.0e-6, nturns=1)
    numpy.testing.assert_equal(pyout, cout)

    set_tilt(pylattice, [1.0e-3], relative=False)
    set_tilt(clattice, [1.0e-3], relative=False)
    pyout = lattice_pass(pylattice, numpy.zeros(6)+1.0e-6, nturns=1)
    cout = lattice_pass(clattice, numpy.zeros(6)+1.0e-6, nturns=1)
    numpy.testing.assert_equal(pyout, cout)


def test_pyintegrator(hmba_lattice):
    params = {'Length': 0,
              'PassMethod': 'pyIdentityPass',
              }
    id_elem = Element('py_id', **params)
    pin = numpy.zeros(6)+1.0e-6
    pout1 = lattice_pass(hmba_lattice, pin.copy(), nturns=1)
    hmba_lattice = hmba_lattice + [id_elem]
    pout2 = lattice_pass(hmba_lattice, pin.copy(), nturns=1)
    numpy.testing.assert_equal(pout1, pout2)
