from __future__ import print_function
import numpy
# noinspection PyUnresolvedReferences
from .atpass import atpass
from .lattice import uint32_refpts


DIMENSION_ERROR = 'Input to lattice_pass() must be a 6xN array.'


def lattice_pass(lattice, r_in, nturns=1, refpts=None, keep_lattice=False):
    """lattice_pass tracks particles through each element of the sequence lattice
    calling the element-specific tracking function specified in the
    lattice[i].PassMethod field.

    Note:

     * lattice_pass(lattice, r_in, refpts=len(line)) is the same as
       lattice_pass(lattice, r_in) since the reference point len(line) is the
       exit of the last element
     * linepass(lattice, r_in, refpts=0) is a copy of r_in since the
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
        6xAxBxC array containing output coordinates of A particles at B
        selected indices for C turns.
    """
    assert r_in.shape[0] == 6 and r_in.ndim in (1, 2), DIMENSION_ERROR
    nparticles = 1 if r_in.ndim == 1 else r_in.shape[1]
    r_in = numpy.asfortranarray(r_in)
    if refpts is None:
        refpts = len(lattice)
    refs = uint32_refpts(refpts, len(lattice))
    result = atpass(lattice, r_in, nturns, refs, int(keep_lattice))
    # atpass returns 6xN array where n = x*y*z;
    # * x is number of particles;
    # * y is number of refpts
    # * z is the number of turns
    # The sequence of output is coordinates for each particle for each refpt
    # for each turn - that is, the first x columns are the x particles at the
    # first refpt on the first # turn, and the first x * y columns are the x
    # particles at all refpts on the first turn.
    # Fortran-order reshaping gathers the elements in this order - from first
    # index of the 4D array to last.
    result = result.reshape((6, nparticles, len(refs), nturns), order='F')
    return result
