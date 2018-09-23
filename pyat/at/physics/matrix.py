"""
transfer matrix related functions

A collection of functions to compute 4x4 and 6x6 transfer matrices
"""

import numpy
from .orbit import *
from ..track import lattice_pass

__all__ = ['find_m44', 'find_m66']

XYDEFSTEP = 6.055454452393343e-006  # Optimal delta?
DPSTEP = 6.055454452393343e-006  # Optimal delta?

# Prepare symplectic identity matrix
_s0 = numpy.zeros((2, 2), order='F')
_s2 = numpy.array([[0, 1], [-1, 0]], order='F')
_s4 = numpy.concatenate((numpy.concatenate((_s2, _s0), axis=1), numpy.concatenate((_s0, _s2), axis=1)), axis=0)
# prepare Identity matrix


def find_m44(ring, dp=0.0, refpts=None, orbit=None, keep_lattice=False, output_orbit=False, **kwargs):
    """find_m44 numerically finds the 4x4 transfer matrix of an accelerator lattice
    for a particle with relative momentum deviation DP

    IMPORTANT!!! find_m44 assumes constant momentum deviation.
    PassMethod used for any element in the lattice SHOULD NOT
    1.  change the longitudinal momentum dP
        (cavities , magnets with radiation, ...)
    2.  have any time dependence (localized impedance, fast kickers, ...)

    m44 = find_m44(lattice, dp=0.0)
        returns a full one-turn matrix at the entrance of the first element
        !!! With this syntax find_m44 assumes that the lattice
        is a ring and first finds the closed orbit

    m44, t = find_m44(lattice, dp, refpts)
        returns 4x4 transfer matrices between the entrance of the first element and each element indexed by refpts.
        t is a 4x4xNrefs array

    m44, t = find_m44(lattice, dp, refpts, full=True)
        Same as above except that matrixes returned in t are full 1-turn matrices at the entrance of each
        element indexed by refpts.

    ... = find_m44(lattice, ..., orbit=closed_orbit)
        Does not search for the closed orbit. Instead closed_orbit,a vector of initial conditions is used.
        This syntax is useful to specify the entrance orbit if lattice is not a ring or to avoid recomputing the
        closed orbit if it is already known.

    m44, t, orbit = find_m44(lattice, ..., output_orbit=True)
        Returns in addition the closed orbit at the entrance of the 1st element

    See also find_m66, find_orbit4
    """
    def mrotate(m):
        m = numpy.squeeze(m)
        return numpy.linalg.multi_dot([m, m44, _s4.T, m.T, _s4])

    xy_step = kwargs.pop('XYStep', XYDEFSTEP)
    full = kwargs.pop('full', False)
    if orbit is None:
        orbit = find_orbit4(ring, dp)
        keep_lattice = True
    # Construct matrix of plus and minus deltas
    dg = numpy.asfortranarray(0.5 * numpy.diag([xy_step] * 6)[:, :4])
    dmat = numpy.concatenate((dg, -dg, numpy.zeros((6, 1))), axis=1)
    # Add the deltas to multiple copies of the closed orbit
    in_mat = orbit.reshape(6, 1) + dmat

    refs = () if refpts is None else refpts
    out_mat = numpy.squeeze(lattice_pass(ring, in_mat, refpts=refs, keep_lattice=keep_lattice), axis=3)
    # out_mat: 8 particles at n refpts for one turn
    # (x + d) - (x - d) / d
    m44 = (in_mat[:4, :4] - in_mat[:4, 4:-1]) / xy_step

    if refpts is not None:
        mstack = (out_mat[:4, :4, :] - out_mat[:4, 4:-1, :]) / xy_step
        if full:
            mstack = numpy.stack(map(mrotate, numpy.split(mstack, mstack.shape[2], axis=2)), axis=2)
        if output_orbit:
            return m44, mstack, out_mat[:, -1, :]
        else:
            return m44, mstack
    else:
        return m44


def find_m66(ring, refpts=None, orbit=None, keep_lattice=False, output_orbit=False, **kwargs):
    """find_m66 numerically finds the 6x6 transfer matrix of an accelerator lattice
    by differentiation of lattice_pass near the closed orbit.
    FINDM66 uses find_orbit6 to search for the closed orbit in 6-D
    In order for this to work the ring MUST have a CAVITY element

    m66 = find_m66(lattice)
        returns the full one-turn 6-by-6 matrix at the entrance of the first element

    m66, t = find_m66(lattice, refpts)
        returns 6x6 transfer matrices between the entrance of the first element and each element indexed by refpts.
        t is 6x6xNrefs array.

    ... = find_m66(lattice, ..., orbit=closed_orbit) - Does not search for closed orbit.
        Does not search for the closed orbit. Instead closed_orbit,a vector of initial conditions is used.
        This syntax is useful to specify the entrance orbit if lattice is not a ring or to avoid recomputing the
        closed orbit if it is already known.

    m66, t, orbit = find_m66(lattice, ..., output_orbit=True)
        Returns in addition the closed orbit at the entrance of the 1st element

    See also find_m44, find_orbit6

    """
    xy_step = kwargs.pop('XYStep', XYDEFSTEP)
    dp_step = kwargs.pop('DPStep', DPSTEP)
    if orbit is None:
        orbit = find_orbit6(ring)
        keep_lattice = True

    # Construct matrix of plus and minus deltas
    scaling = numpy.array([xy_step, xy_step, xy_step, xy_step, dp_step, dp_step])
    dg = numpy.asfortranarray(0.5*numpy.diag(scaling))
    dmat = numpy.concatenate((dg, -dg, numpy.zeros((6, 1))), axis=1)

    in_mat = orbit.reshape(6, 1) + dmat

    refs = () if refpts is None else refpts
    out_mat = numpy.squeeze(lattice_pass(ring, in_mat, refpts=refs, keep_lattice=keep_lattice), axis=3)
    # out_mat: 12 particles at n refpts for one turn
    # (x + d) - (x - d) / d
    m66 = (in_mat[:, :6] - in_mat[:, 6:-1]) / scaling.reshape((1, 6))

    if refpts is not None:
        mstack = (out_mat[:, :6, :] - out_mat[:, 6:-1, :]) / xy_step
        if output_orbit:
            return m66, mstack, out_mat[:, -1, :]
        else:
            return m66, mstack
    else:
        return m66
