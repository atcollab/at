"""
transfer matrix related functions

A collection of functions to compute 4x4 and 6x6 transfer matrices
"""

import numpy
import at
import math
from .orbit import *

__all__ = ['find_m44', 'find_m66', 'get_twiss']

XYDEFSTEP = 6.055454452393343e-006  # Optimal delta?
DPSTEP = 6.055454452393343e-006  # Optimal delta?
DDP = 1e-8

# dtype for structured array containing Twiss parameters - see get_twiss()
TWISS_DTYPE = [('idx', numpy.uint32),
               ('s_pos', numpy.float64),
               ('closed_orbit', numpy.float64, (6,)),
               ('dispersion', numpy.float64, (4,)),
               ('alpha', numpy.float64, (2,)),
               ('beta', numpy.float64, (2,)),
               ('mu', numpy.float64, (2,)),
               ('m44', numpy.float64, (4, 4))]

# Prepare symplectic identity matrix
_s0 = numpy.zeros((2, 2), order='F')
_s2 = numpy.array([[0, 1], [-1, 0]], order='F')
_s4 = numpy.concatenate((numpy.concatenate((_s2, _s0), axis=1), numpy.concatenate((_s0, _s2), axis=1)), axis=0)
# prepare Identity matrix


def find_m44(ring, dp=0.0, refpts=None, orbit=None, output_orbit=False, **kwargs):
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
    keeplattice = False
    if orbit is None:
        orbit = find_orbit4(ring, dp)
        keeplattice = True
    # Construct matrix of plus and minus deltas
    dg = numpy.asfortranarray(0.5 * numpy.diag([xy_step] * 6)[:, :4])
    dmat = numpy.concatenate((dg, -dg, numpy.zeros((6, 1))), axis=1)
    # Add the deltas to multiple copies of the closed orbit
    in_mat = orbit.reshape(6, 1) + dmat

    refs = () if refpts is None else refpts
    out_mat = numpy.squeeze(at.lattice_pass(ring, in_mat, refpts=refs, keep_lattice=keeplattice), axis=3)
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


def find_m66(ring, refpts=None, orbit=None, output_orbit=False, **kwargs):
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
    keeplattice = False
    if orbit is None:
        orbit = find_orbit6(ring)
        keeplattice = True

    # Construct matrix of plus and minus deltas
    scaling = numpy.array([xy_step, xy_step, xy_step, xy_step, dp_step, dp_step])
    dg = numpy.asfortranarray(0.5*numpy.diag(scaling))
    dmat = numpy.concatenate((dg, -dg, numpy.zeros((6, 1))), axis=1)

    in_mat = orbit.reshape(6, 1) + dmat

    refs = () if refpts is None else refpts
    out_mat = numpy.squeeze(at.lattice_pass(ring, in_mat, refpts=refs, keep_lattice=keeplattice), axis=3)
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


def betatron_phase_unwrap(m):
    """
    Unwrap negative jumps in betatron.
    """
    dp = numpy.diff(m)
    jumps = numpy.append([0], dp) < 0
    return m + numpy.cumsum(jumps) * numpy.pi


def get_twiss(ring, dp=0.0, refpts=None, get_chrom=False, ddp=DDP):
    """Determine Twiss parameters by first finding the transfer matrix.

    The Twiss structured array has nrefs elements, so:
     * twiss['idx'].shape is (nrefs,)
     * twiss['closed_orbit'].shape is (nrefs, 6).

    Returns:
        twiss - structured array
        tune - numpy array of shape (2,)
        chrom - numpy array of shape (2,)
    """
    def twiss22(mat, ms):
        """
        Calculate Twiss parameters from the standard 2x2 transfer matrix
        (i.e. x or y).
        """
        sin_mu_end = (numpy.sign(mat[0, 1]) *
                      math.sqrt(-mat[0, 1] * mat[1, 0] -
                                (mat[0, 0] - mat[1, 1]) ** 2 / 4))
        alpha0 = (mat[0, 0] - mat[1, 1]) / 2.0 / sin_mu_end
        beta0 = mat[0, 1] / sin_mu_end
        beta = ((ms[0, 0, :] * beta0 - ms[0, 1, :] * alpha0) **
                2 + ms[0, 1, :] ** 2) / beta0
        alpha = -((ms[0, 0, :] * beta0 - ms[0, 1, :] * alpha0) *
                  (ms[1, 0, :] * beta0 - ms[1, 1, :] * alpha0) +
                  ms[0, 1, :] * ms[1, 1, :]) / beta0
        mu = numpy.arctan(ms[0, 1, :] / (ms[0, 0, :] * beta0 - ms[0, 1, :] * alpha0))
        mu = betatron_phase_unwrap(mu)
        return alpha, beta, mu

    chrom = None

    refpts = at.uint32_refpts(refpts, len(ring))
    nrefs = refpts.size
    if refpts[-1] != len(ring):
        refpts = numpy.append(refpts, [len(ring)])

    orbit4, orbit = find_orbit4(ring, dp, refpts)
    m44, mstack = find_m44(ring, dp, refpts, orbit=orbit4)

    ax, bx, mx = twiss22(m44[:2, :2], mstack[:2, :2, :])
    ay, by, my = twiss22(m44[2:, 2:], mstack[2:, 2:, :])

    tune = numpy.array((mx[-1], my[-1])) / (2 * numpy.pi)
    twiss = numpy.zeros(nrefs, dtype=TWISS_DTYPE)
    twiss['idx'] = refpts[:nrefs]
    # Use rollaxis to get the arrays in the correct shape for the twiss
    # structured array - that is, with nrefs as the first dimension.
    twiss['s_pos'] = at.get_s_pos(ring, refpts[:nrefs])
    twiss['closed_orbit'] = numpy.rollaxis(orbit, -1)[:nrefs]
    twiss['m44'] = numpy.rollaxis(mstack, -1)[:nrefs]
    twiss['alpha'] = numpy.rollaxis(numpy.vstack((ax, ay)), -1)[:nrefs]
    twiss['beta'] = numpy.rollaxis(numpy.vstack((bx, by)), -1)[:nrefs]
    twiss['mu'] = numpy.rollaxis(numpy.vstack((mx, my)), -1)[:nrefs]
    twiss['dispersion'] = numpy.NaN
    # Calculate chromaticity by calling this function again at a slightly
    # different momentum.
    if get_chrom:
        twissb, tuneb, _ = get_twiss(ring, dp + ddp, refpts[:nrefs])
        chrom = (tuneb - tune) / ddp
        twiss['dispersion'] = (twissb['closed_orbit'] - twiss['closed_orbit'])[:, :4] / ddp

    return twiss, tune, chrom
