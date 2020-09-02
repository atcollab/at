"""Simple tests to make sure integrators load in Python."""
import numpy
import pytest

from at import atpass, elements
from at.lattice import Lattice


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
