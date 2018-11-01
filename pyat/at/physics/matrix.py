"""
transfer matrix related functions

A collection of functions to compute 4x4 and 6x6 transfer matrices
"""

import numpy
from scipy.linalg import block_diag
from ..physics import find_orbit4, find_orbit6
from ..tracking import lattice_pass, element_pass

__all__ = ['find_m44', 'find_m66', 'find_elem_m66', 'jmat']

XYDEFSTEP = 6.055454452393343e-006  # Optimal delta?
DPSTEP = 6.055454452393343e-006  # Optimal delta?

# Prepare symplectic identity matrix
_j2 = numpy.array([[0., 1.], [-1., 0.]])
_jm = [_j2, block_diag(_j2, _j2), block_diag(_j2, _j2, _j2)]


def jmat(ind):
    """
    Return the antisymetric block diagonal matrix [[0, 1][-1, 0]]

    INPUT
        ind     1, 2 or 3. Matrix dimension

    OUTPUT
        jm      block diagonal matrix, (2, 2) or (4, 4) or (6, 6)
    """
    return _jm[ind - 1]


_jmt = jmat(2)


def find_m44(ring, dp=0.0, refpts=None, orbit=None, keep_lattice=False, **kwargs):
    """find_m44 numerically finds the 4x4 transfer matrix of an accelerator lattice
    for a particle with relative momentum deviation DP

    IMPORTANT!!! find_m44 assumes constant momentum deviation.
    PassMethod used for any element in the lattice SHOULD NOT
    1.  change the longitudinal momentum dP
        (cavities , magnets with radiation, ...)
    2.  have any time dependence (localized impedance, fast kickers, ...)

    m44, t = find_m44(lattice, dp=0.0, refpts)
        return 4x4 transfer matrices between the entrance of the first element and each element indexed by refpts.
        m44:    full one-turn matrix at the entrance of the first element
        t:      4x4 transfer matrices between the entrance of the first element and each element indexed by refpts:
                (Nrefs, 4, 4) array

        Unless an input orbit is introduced, find_m44 assumes that the lattice is a ring and first
        finds the closed orbit.

    KEYWORDS
        keep_lattice=False  When True, assume no lattice change since the previous tracking.
        full=False          When True, matrices are full 1-turn matrices at the entrance of each
                            element indexed by refpts.
        orbit=None          Avoids looking for the closed orbit if is already known ((6,) array)
        XYStep=6.055e-6     transverse step for numerical computation

    See also find_m66, find_orbit4
    """

    def mrotate(m):
        m = numpy.squeeze(m)
        return numpy.linalg.multi_dot([m, m44, _jmt.T, m.T, _jmt])

    xy_step = kwargs.pop('XYStep', XYDEFSTEP)
    full = kwargs.pop('full', False)
    if orbit is None:
        orbit, _ = find_orbit4(ring, dp, keep_lattice=keep_lattice)
        keep_lattice = True
    # Construct matrix of plus and minus deltas
    dg = numpy.asfortranarray(0.5 * numpy.diag([xy_step] * 6)[:, :4])
    dmat = numpy.concatenate((dg, -dg), axis=1)
    # Add the deltas to multiple copies of the closed orbit
    in_mat = orbit.reshape(6, 1) + dmat

    refs = () if refpts is None else refpts
    out_mat = numpy.rollaxis(
              numpy.squeeze(lattice_pass(ring, in_mat, refpts=refs, keep_lattice=keep_lattice), axis=3), -1)
    # out_mat: 8 particles at n refpts for one turn
    # (x + d) - (x - d) / d
    m44 = (in_mat[:4, :4] - in_mat[:4, 4:]) / xy_step

    if len(refs) > 0:
        mstack = (out_mat[:, :4, :4] - out_mat[:, :4, 4:]) / xy_step
        if full:
            # mstack = numpy.stack(map(mrotate, numpy.split(mstack, mstack.shape[0], axis=0)), axis=0)
            mstack = numpy.stack((mrotate(mat) for mat in mstack), axis=0)
    else:
        mstack = numpy.empty((0, 4, 4), dtype=float)

    return m44, mstack


def find_m66(ring, refpts=None, orbit=None, keep_lattice=False, **kwargs):
    """find_m66 numerically finds the 6x6 transfer matrix of an accelerator lattice
    by differentiation of lattice_pass near the closed orbit.
    find_m66 uses find_orbit6 to search for the closed orbit in 6-D
    In order for this to work the ring MUST have a CAVITY element

    m66, t = find_m66(lattice, refpts)
        m66:    full one-turn 6-by-6 matrix at the entrance of the first element
        t:      6x6 transfer matrices between the entrance of the first element and each element indexed by refpts:
                (nrefs, 6, 6) array.

    KEYWORDS
        keep_lattice=False  When True, assume no lattice change since the previous tracking.
        orbit=None          Avoids looking for the closed orbit if is already known ((6,) array)
        XYStep=6.055e-6     transverse step for numerical computation
        DPStep=6.055e-6     longitudinal step for numerical computation

    See also find_m44, find_orbit6

    """
    xy_step = kwargs.pop('XYStep', XYDEFSTEP)
    dp_step = kwargs.pop('DPStep', DPSTEP)
    if orbit is None:
        orbit, _ = find_orbit6(ring, keep_lattice=keep_lattice)
        keep_lattice = True

    # Construct matrix of plus and minus deltas
    scaling = numpy.array([xy_step, xy_step, xy_step, xy_step, dp_step, dp_step])
    dg = numpy.asfortranarray(0.5 * numpy.diag(scaling))
    dmat = numpy.concatenate((dg, -dg), axis=1)

    in_mat = orbit.reshape(6, 1) + dmat

    refs = () if refpts is None else refpts
    out_mat = numpy.rollaxis(
              numpy.squeeze(lattice_pass(ring, in_mat, refpts=refs, keep_lattice=keep_lattice), axis=3), -1)
    # out_mat: 12 particles at n refpts for one turn
    # (x + d) - (x - d) / d
    m66 = (in_mat[:, :6] - in_mat[:, 6:]) / scaling.reshape((1, 6))

    if len(refs) > 0:
        mstack = (out_mat[:, :, :6] - out_mat[:, :, 6:]) / xy_step
    else:
        mstack = numpy.empty((0, 6, 6), dtype=float)

    return m66, mstack


def find_elem_m66(elem, orbit=None, **kwargs):
    """
    Numerically find the 6x6 transfer matrix of a single element

    INPUT
        elem                AT element

    KEYWORDS
        orbit=None          closed orbit at the entrance of the element. Default: 0.0
        XYStep=6.055e-6     transverse step for numerical computation
        DPStep=6.055e-6     longitudinal step for numerical computation

    OUTPUT
        m66                 (6, 6) transfer matrix
    """
    xy_step = kwargs.pop('XYStep', XYDEFSTEP)
    dp_step = kwargs.pop('DPStep', DPSTEP)
    if orbit is None:
        orbit = numpy.zeros((6,))

    # Construct matrix of plus and minus deltas
    scaling = numpy.array([xy_step, xy_step, xy_step, xy_step, dp_step, dp_step])
    dg = numpy.asfortranarray(0.5 * numpy.diag(scaling))
    dmat = numpy.concatenate((dg, -dg), axis=1)

    in_mat = orbit.reshape(6, 1) + dmat
    element_pass(elem, in_mat)
    m66 = (in_mat[:, :6] - in_mat[:, 6:]) / scaling.reshape((1, 6))
    return m66
