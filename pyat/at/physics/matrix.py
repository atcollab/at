"""
transfer matrix related functions

A collection of functions to compute 4x4 and 6x6 transfer matrices
"""
import numpy
from ..lattice import Lattice, Element, DConstant, Refpts, Orbit
from ..lattice import frequency_control, get_uint32_index
from ..lattice.elements import Dipole, M66
from ..tracking import internal_lpass, internal_epass
from .orbit import find_orbit4, find_orbit6
from .amat import jmat, symplectify

__all__ = ['find_m44', 'find_m66', 'find_elem_m66', 'gen_m66_elem']

_jmt = jmat(2)


def find_m44(ring: Lattice, dp: float = None, refpts: Refpts = None,
             dct: float = None, df: float = None,
             orbit: Orbit = None, keep_lattice: bool = False, **kwargs):
    """One turn 4x4 transfer matrix

    :py:func:`find_m44` finds the 4x4 transfer matrix of an accelerator
    lattice by differentiation of trajectories near the closed orbit.

    Important:
        :py:func:`find_m44` assumes constant momentum deviation.
        The ``PassMethod`` used for any element **SHOULD NOT**:

        1.  change the longitudinal momentum dP
            (cavities , magnets with radiation, ...)
        2.  have any time dependence (localized impedance, fast kickers, ...)

    Unless an input orbit is introduced, :py:func:`find_m44` assumes that the
    lattice is a ring and first finds the closed orbit.

    Parameters:
        ring:           Lattice description (radiation must be OFF)
        dp:             Momentum deviation.
        refpts:         Observation points
        dct:            Path lengthening.
        df:             Deviation of RF frequency.
        orbit:          Avoids looking for initial the closed orbit if it is
          already known ((6,) array).
        keep_lattice:   Assume no lattice change since the previous tracking.
          Default: :py:obj:`False`

    Keyword Args:
        full (bool):    If :py:obj:`True`, matrices are full 1-turn matrices
          at the entrance of each element indexed by refpts. If :py:obj:`False`
          (Default), matrices are between the entrance of the first element and
          the entrance of the selected element
        XYStep (float): Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`

    Returns:
        m44:    full one-turn matrix at the entrance of the first element
        ms:     4x4 transfer matrices between the entrance of the first
          element and each element indexed by refpts: (Nrefs, 4, 4) array

    See also:
         :py:func:`find_m66`, :py:func:`.find_orbit4`
    """

    def mrotate(m):
        m = numpy.squeeze(m)
        return m.dot(m44.dot(_jmt.T.dot(m.T.dot(_jmt))))

    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    full = kwargs.pop('full', False)
    if orbit is None:
        orbit, _ = find_orbit4(ring, dp=dp, dct=dct, df=df,
                               keep_lattice=keep_lattice, XYStep=xy_step)
        keep_lattice = True
    # Construct matrix of plus and minus deltas
    # scaling = 2*xy_step*numpy.array([1.0, 0.1, 1.0, 0.1])
    scaling = xy_step * numpy.array([1.0, 1.0, 1.0, 1.0])
    dg = numpy.asfortranarray(
        numpy.concatenate((0.5 * numpy.diag(scaling), numpy.zeros((2, 4)))))
    dmat = numpy.concatenate((dg, -dg), axis=1)
    # Add the deltas to multiple copies of the closed orbit
    in_mat = orbit.reshape(6, 1) + dmat

    refs = get_uint32_index(ring, refpts)
    out_mat = numpy.rollaxis(
        numpy.squeeze(internal_lpass(ring, in_mat, refpts=refs,
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


@frequency_control
def find_m66(ring: Lattice, refpts: Refpts = None,
             orbit: Orbit = None, keep_lattice: bool = False, **kwargs):
    """One-turn 6x6 transfer matrix

    :py:func:`find_m66` finds the 6x6 transfer matrix of an accelerator
    lattice by differentiation of trajectories near the closed orbit.

    Parameters:
        ring:           Lattice description
        refpts:         Observation points
        orbit:          Avoids looking for initial the closed orbit if it is
          already known ((6,) array).
        keep_lattice:   Assume no lattice change since the previous tracking.
          Default: :py:obj:`False`

    Keyword Args:
        XYStep (float)  Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float): Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`

    Returns:
        m66:    full one-turn matrix at the entrance of the first element
        ms:     6x6 transfer matrices between the entrance of the first
          element and each element indexed by refpts: (Nrefs, 6, 6) array

    See also:
         :py:func:`find_m44`, :py:func:`.find_orbit6`
    """
    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    dp_step = kwargs.pop('DPStep', DConstant.DPStep)
    if orbit is None:
        if ring.radiation:
            orbit, _ = find_orbit6(ring, keep_lattice=keep_lattice,
                                   XYStep=xy_step, DPStep=dp_step, **kwargs)
        else:
            orbit, _ = find_orbit4(ring, keep_lattice=keep_lattice,
                                   XYStep=xy_step, **kwargs)
        keep_lattice = True

    # Construct matrix of plus and minus deltas
    # scaling = 2*xy_step*numpy.array([1.0, 0.1, 1.0, 0.1, 1.0, 1.0])
    scaling = xy_step * numpy.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0]) + \
        dp_step * numpy.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
    dg = numpy.asfortranarray(0.5 * numpy.diag(scaling))
    dmat = numpy.concatenate((dg, -dg), axis=1)

    in_mat = orbit.reshape(6, 1) + dmat

    refs = get_uint32_index(ring, refpts)
    out_mat = numpy.rollaxis(
        numpy.squeeze(internal_lpass(ring, in_mat, refpts=refs,
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


def find_elem_m66(elem: Element, orbit: Orbit = None, **kwargs):
    """Single element 6x6 transfer matrix

    Numerically finds the 6x6 transfer matrix of a single element

    Parameters:
        elem:           AT element
        orbit:          Closed orbit at the entrance of the element,
          default: 0.0

    Keyword Args:
        XYStep (float): Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        particle (Particle):    circulating particle.
          Default: :code:`lattice.particle` if existing,
          otherwise :code:`Particle('relativistic')`
        energy (float):         lattice energy. Default 0.

    Returns:
        m66:            6x6 transfer matrix
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
    internal_epass(elem, in_mat, **kwargs)
    m66 = (in_mat[:, :6] - in_mat[:, 6:]) / scaling
    return m66


def gen_m66_elem(ring: Lattice, 
                 ringrad: Lattice = None,
                 o6b: Orbit = None,
                 o6e: Orbit = None,
                 o6brad: Orbit = None,
                 o6erad: Orbit = None,
                 **kwargs) -> M66:
    """Converts a ring to a linear 6x6 matrix. This element handles 2 passmethods
    ``Matrix66Pass`` and ``Matrix66RadPass``, that can be activated within a lattice
    using :py:func:`enable_6()` or :py:func:`disable_6d`.
    This provides flexibility to turn ON and OFF the 6D motion or simply the radiation
    damping, depending on the configuration of the rings provided as input.

    Parameters:
        ring:       Lattice description. Is used by ``Matrix66Pass``
        ringrad:    Optional lattice with radiations. Is used by ``Matrix66RadPass``
        o6b:        entrace orbit without radiations, use by ``Matrix66Pass``
        o6e:        exit orbit without radiations, use by ``Matrix66Pass``
        o6b:        entrace orbit with radiations, use by ``Matrix66RadPass``
        o6e:        exit orbit with radiations, use by ``Matrix66RadPass``

    Returns:
        m66:        :py:obj:`M66` object
    """
    length = ring.circumference
    kwargs.update({"Length": length})
    m66_mat, _ = find_m66(ring, [], orbit=o6b)
    m66_mat_rad, _ = find_m66(ringrad, [], orbit=o6brad)
    kwargs.update({"M66Rad": m66_mat_rad})
    if o6b is not None:
        kwargs.update({"T1": -o6b})
    if o6e is not None:
        kwargs.update({"T2": o6e})
    if o6brad is not None:
        kwargs.update({"T1rad": -o6brad})
    if o6erad is not None:
        kwargs.update({"T2rad": o6erad})  
    lin_elem = M66('Linear', m66_mat, **kwargs)
    return lin_elem
       

Lattice.find_m44 = find_m44
Lattice.find_m66 = find_m66
