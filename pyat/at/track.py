from __future__ import print_function
import numpy
# noinspection PyUnresolvedReferences
from .atpass import atpass
from .lattice import uint32_refpts


def lattice_pass(lattice, r_in, nturns=1, refpts=None, keep_lattice=False):
    """pass tracks particles through each element of the sequence lattice
    calling the element-specific tracking function specified in the
    lattice[i].PassMethod field.

    Note:

     * lattice_pass(lattice, r_in, refpts=len(line)) is the same as
       lattice_pass(lattice, r_in) since the reference point len(line) is the
       exit of the last element
     * linepass(lattice, r_in, refpts=0]) is a copy of r_in since the
       reference point 0 is the entrance of the first element

    Args:
        lattice: sequence of AT elements
        r_in: 6xN array: input coordinates of N particles
        nturns: number of passes through the lattice line
        refpts: indices of elements at which to return coordinates (see
                lattice.py)
        keep_lattice: use elements persisted from a previous call to at.atpass.
                      If True, assume that the lattice has not changed since
                      that previous call.

    Returns:
        6xN array containing output coordinates of x particles at y selected
        indices for z turns; N = x * y * z. The sequence of output is
        coordinates for each particle for each refpt for each turn - that is,
        the first x columns are the x particles at the first refpt on the first
        turn, and the first x * y columns are the x particles at all refpts on
        the first turn.
    """
    assert r_in.shape[0] == 6
    r_in = numpy.asfortranarray(r_in)
    if refpts is None:
        refpts = len(lattice)
    refs = uint32_refpts(refpts, len(lattice))
    return atpass(lattice, r_in, nturns, refs, int(keep_lattice))
