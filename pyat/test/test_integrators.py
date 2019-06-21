"""Simple tests to make sure integrators load in Python."""
import pytest

from at import atpass, elements
from at.lattice import Lattice


def test_exact_hamiltonian_pass(rin):
    drift = elements.Multipole('m1', 1, [0, 0, 0, 0], [0, 0, 0, 0])
    drift.Type = 0
    drift.PassMethod = 'ExactHamiltonianPass'
    drift.BendingAngle = 0
    l = Lattice([drift], name='lat', energy=5)
    atpass(l, rin, 1)


@pytest.mark.parametrize('passmethod', ('GWigSymplecticPass', 'GWigSymplecticRadPass'))
def test_gwig_symplectic_pass(rin, passmethod):
    # Parameters copied from one of the Diamond wigglers.
    wiggler = elements.Wiggler('w', 1.15, 0.05, 0, 5, 4, [1, 1, 0, 1, 1, 0], [], 3e9)
    wiggler.PassMethod = passmethod
    l = Lattice([wiggler], name='lat', energy=5)
    atpass(l, rin, 1)