import numpy
from at.tracking import atpass, elempass
from at.lattice import uint32_refpts


__all__ = ['lattice_pass', 'element_pass']

DIMENSION_ERROR = 'Input to lattice_pass() must be a 6xN array.'


def lattice_pass(lattice, r_in, nturns=1, refpts=None, keep_lattice=False):
    """lattice_pass tracks particles through each element of a lattice
    calling the element-specific tracking function specified in the
    lattice[i].PassMethod field.

    Note:

     * lattice_pass(lattice, r_in, refpts=len(line)) is the same as
       lattice_pass(lattice, r_in) since the reference point len(line) is the
       exit of the last element
     * lattice_pass(lattice, r_in, refpts=0) is a copy of r_in since the
       reference point 0 is the entrance of the first element

    PARAMETERS
        lattice:    iterable of AT elements
        r_in:       (6, N) array: input coordinates of N particles.
                    r_in is modified in-place and reports the coordinates at
                    the end of the tracking. For the the best efficiency, r_in
                    should be given as F_CONTIGUOUS numpy array.
        nturns:     number of passes through the lattice line
        refpts      elements at which data is returned. It can be:
                    1) an integer in the range [-len(ring), len(ring)-1]
                       selecting the element according to python indexing
                       rules. As a special case, len(ring) is allowed and
                       refers to the end of the last element,
                    2) an ordered list of such integers without duplicates,
                    3) a numpy array of booleans of maximum length
                       len(ring)+1, where selected elements are True.
                    Defaults to None, meaning no refpts, equivelent to
                    passing an empty array for calculation purposes.
        keep_lattice: use elements persisted from a previous call to at.atpass.
                    If True, assume that the lattice has not changed since
                    that previous call.

    OUTPUT
        (6, N, R, T) array containing output coordinates of N particles
        at R reference points for T turns.
    """
    assert r_in.shape[0] == 6 and r_in.ndim in (1, 2), DIMENSION_ERROR
    if not isinstance(lattice, list):
        lattice = list(lattice)
    nelems = len(lattice)
    if refpts is None:
        refpts = nelems
    refs = uint32_refpts(refpts, nelems)
    # atpass returns 6xAxBxC array where n = x*y*z;
    # * A is number of particles;
    # * B is number of refpts
    # * C is the number of turns
    if r_in.flags.f_contiguous:
        return atpass(lattice, r_in, nturns, refs, int(keep_lattice))
    else:
        r_fin = numpy.asfortranarray(r_in)
        r_out = atpass(lattice, r_fin, nturns, refs, int(keep_lattice))
        r_in[:] = r_fin[:]
        return r_out


def element_pass(element, r_in):
    r_in = numpy.asfortranarray(r_in)
    elempass(element, r_in)
