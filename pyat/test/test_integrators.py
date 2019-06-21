"""Simple tests to make sure integrators load in Python."""

from at import atpass, elements
from at.lattice import Lattice


def test_exact_hamiltonian_pass(rin):
    drift = elements.Multipole('m1', 1, [0, 0, 0, 0], [0, 0, 0, 0])
    drift.Type = 0
    drift.PassMethod = 'ExactHamiltonianPass'
    drift.BendingAngle = 0
    l = Lattice([drift], name='lat', energy=5)
    atpass(l, rin, 1)
