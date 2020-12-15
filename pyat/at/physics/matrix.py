"""
transfer matrix related functions

A collection of functions to compute 4x4 and 6x6 transfer matrices
"""

import numpy
from at.lattice import Lattice, uint32_refpts, get_refpts
from at.lattice.elements import Bend, M66
from at.tracking import lattice_pass, element_pass
from at.physics import find_orbit4, find_orbit6, jmat, symplectify, DConstant

__all__ = ['find_m44', 'find_m66', 'find_elem_m66', 'gen_m66_elem']

_jmt = jmat(2)


def find_m44(ring, dp=0.0, refpts=None, orbit=None, keep_lattice=False,
             **kwargs):
    """find_m44 numerically finds the 4x4 transfer matrix of an accelerator
    lattice for a particle with relative momentum deviation DP

    IMPORTANT!!! find_m44 assumes constant momentum deviation.
    PassMethod used for any element in the lattice SHOULD NOT
    1.  change the longitudinal momentum dP
        (cavities , magnets with radiation, ...)
    2.  have any time dependence (localized impedance, fast kickers, ...)

    m44, t = find_m44(lattice, dp=0.0, refpts)
        return 4x4 transfer matrices between the entrance of the first element
        and each element indexed by refpts.
            m44:    full one-turn matrix at the entrance of the first element
            t:      4x4 transfer matrices between the entrance of the first
                    element and each element indexed by refpts:
                    (Nrefs, 4, 4) array

    Unless an input orbit is introduced, find_m44 assumes that the lattice is
    a ring and first finds the closed orbit.

    PARAMETERS
        ring            lattice description
        dp              momentum deviation. Defaults to 0
        refpts          elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.
                        Defaults to None, if refpts is None an empty array is
                        returned for mstack.

    KEYWORDS
        keep_lattice=False  When True, assume no lattice change since the
                            previous tracking.
        full=False          When True, matrices are full 1-turn matrices at
                            the entrance of each
                            element indexed by refpts.
        orbit=None          Avoids looking for the closed orbit if is already
                            known (6,) array
        XYStep=1.e-8        transverse step for numerical computation

    See also find_m66, find_orbit4
    """

    def mrotate(m):
        m = numpy.squeeze(m)
        return m.dot(m44.dot(_jmt.T.dot(m.T.dot(_jmt))))

    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    full = kwargs.pop('full', False)
    if orbit is None:
        orbit, _ = find_orbit4(ring, dp, keep_lattice=keep_lattice,
                               XYStep=xy_step)
        keep_lattice = True
    # Construct matrix of plus and minus deltas
    # scaling = 2*xy_step*numpy.array([1.0, 0.1, 1.0, 0.1])
    scaling = xy_step * numpy.array([1.0, 1.0, 1.0, 1.0])
    dg = numpy.asfortranarray(
        numpy.concatenate((0.5 * numpy.diag(scaling), numpy.zeros((2, 4)))))
    dmat = numpy.concatenate((dg, -dg), axis=1)
    # Add the deltas to multiple copies of the closed orbit
    in_mat = orbit.reshape(6, 1) + dmat

    refs = uint32_refpts(refpts, len(ring))
    out_mat = numpy.rollaxis(
        numpy.squeeze(lattice_pass(ring, in_mat, refpts=refs,
                                   keep_lattice=keep_lattice), axis=3), -1
    )
    # out_mat: 8 particles at n refpts for one turn
    # (x + d) - (x - d) / d
    m44 = (in_mat[:4, :4] - in_mat[:4, 4:]) / scaling

    if len(refs) > 0:
        mstack = (out_mat[:, :4, :4] - out_mat[:, :4, 4:]) / scaling
        if full:
            mstack = numpy.stack([mrotate(mat) for mat in mstack], axis=0)
    else:
        mstack = numpy.empty((0, 4, 4), dtype=float)

    return m44, mstack


def find_m66(ring, refpts=None, orbit=None, keep_lattice=False, **kwargs):
    """find_m66 numerically finds the 6x6 transfer matrix of an accelerator
    lattice by differentiation of lattice_pass near the closed orbit.
    find_m66 uses find_orbit6 to search for the closed orbit in 6-D
    In order for this to work the ring MUST have a CAVITY element

    m66, t = find_m66(lattice, refpts)
        m66:    full one-turn 6-by-6 matrix at the entrance of the
                first element.
        t:      6x6 transfer matrices between the entrance of the first
                element and each element indexed by refpts (nrefs, 6, 6) array.

    PARAMETERS
        ring            lattice description
        dp              momentum deviation. Defaults to 0
        refpts          elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.
                        Defaults to None, if refpts is None an empty array is
                        returned for mstack.

    KEYWORDS
        keep_lattice=False  When True, assume no lattice change since the
                            previous tracking.
        orbit=None          Avoids looking for the closed orbit if is already
                            known (6,) array
        XYStep=1.e-8        transverse step for numerical computation

    See also find_m44, find_orbit6
    """
    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    dp_step = kwargs.pop('DPStep', DConstant.DPStep)
    if orbit is None:
        orbit, _ = find_orbit6(ring, keep_lattice=keep_lattice,
                               XYStep=xy_step, DPStep=dp_step)
        keep_lattice = True

    # Construct matrix of plus and minus deltas
    # scaling = 2*xy_step*numpy.array([1.0, 0.1, 1.0, 0.1, 1.0, 1.0])
    scaling = xy_step * numpy.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0]) + \
              dp_step * numpy.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
    dg = numpy.asfortranarray(0.5 * numpy.diag(scaling))
    dmat = numpy.concatenate((dg, -dg), axis=1)

    in_mat = orbit.reshape(6, 1) + dmat

    refs = uint32_refpts(refpts, len(ring))
    out_mat = numpy.rollaxis(
        numpy.squeeze(lattice_pass(ring, in_mat, refpts=refs,
                                   keep_lattice=keep_lattice), axis=3), -1
    )
    # out_mat: 12 particles at n refpts for one turn
    # (x + d) - (x - d) / d
    m66 = (in_mat[:, :6] - in_mat[:, 6:]) / scaling

    if len(refs) > 0:
        mstack = (out_mat[:, :, :6] - out_mat[:, :, 6:]) / scaling
    else:
        mstack = numpy.empty((0, 6, 6), dtype=float)

    return m66, mstack


def find_elem_m66(elem, orbit=None, **kwargs):
    """
    Numerically find the 6x6 transfer matrix of a single element

    INPUT
        elem                AT element

    KEYWORDS
        orbit=None          closed orbit at the entrance of the element,
                            default: 0.0
        XYStep=1.e-8        transverse step for numerical computation

    OUTPUT
        m66                 (6, 6) transfer matrix
    """
    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    if orbit is None:
        orbit = numpy.zeros((6,))

    # Construct matrix of plus and minus deltas
    # scaling = 2*xy_step*numpy.array([1.0, 0.1, 1.0, 0.1, 1.0, 1.0])
    scaling = xy_step * numpy.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    dg = numpy.asfortranarray(0.5 * numpy.diag(scaling))
    dmat = numpy.concatenate((dg, -dg), axis=1)

    in_mat = orbit.reshape(6, 1) + dmat
    element_pass(elem, in_mat)
    m66 = (in_mat[:, :6] - in_mat[:, 6:]) / scaling
    return m66


def gen_m66_elem(ring, o4b, o4e):
    """
    converts a ring to a linear 6x6 matrix tracking elemtn
    """

    dip_inds = get_refpts(ring, Bend)
    theta = numpy.array([ring[ind].BendingAngle for ind in dip_inds])
    lendp = numpy.array([ring[ind].Length for ind in dip_inds])
    s_pos = ring.get_s_pos()
    s = numpy.diff(numpy.array([s_pos[0], s_pos[-1]]))[0]
    i2 = numpy.sum(numpy.abs(theta * theta / lendp))
    m66_mat, _ = find_m66(ring, [], o4b)
    if ring.radiation is False:
        m66_mat = symplectify(m66_mat)  # remove for damping
    lin_elem = M66('Linear', m66_mat, T1=-o4b, T2=o4e, Length=s, I2=i2)
    return lin_elem


Lattice.find_m44 = find_m44
Lattice.find_m66 = find_m66
