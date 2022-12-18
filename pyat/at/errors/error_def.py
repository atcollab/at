import numpy as np
from typing import Optional, Union
from at.lattice import refpts_iterator, Lattice, Refpts
from at.lattice import shift_elem, rotate_elem
from at.physics import find_orbit, get_optics
from at.tracking import lattice_pass
from scipy.stats import truncnorm, norm


__all__ = ['find_orbit_err', 'get_optics_err', 'lattice_pass_err',
           'assign_errors', 'enable_errors']

_BPM_ATTRS = ('BPMGain','BPMOffset','BPMTilt')
_ERR_ATTRS = ('PolynomBErr', 'PolynomAErr', 'ShiftErr', 'RotationErr')


def _truncated_randn(truncation=None, **kwargs):
    if truncation is not None:
        return truncnorm.rvs(-truncation, truncation, **kwargs)
    else:
        return norm.rvs(**kwargs)

def assign_errors(ring: Lattice, key, truncation: Optional[float] = None,
                  seed: Optional[Union[int, np.random.Generator]] = None,
                  **kwargs):
    r"""Assign errors to selected elements

    The errors are stored as additional attributes but are not enabled for
    tracking (see :py:func:`enable_errors`).

    Args:
        ring:   Lattice description.
        key:    Element selector. May be:

                  #. an element instance, will return all elements of the same
                     type in the lattice, e.g. :pycode:`Drift('d1', 1.0)`
                  #. an element type, will return all elements of that type in
                     the lattice, e.g. :pycode:`at.Sextupole`
                  #. a string to match against elements' ``FamName``, supports
                     Unix shell-style wildcards, e.g. ``'BPM_*1'``
        truncation:     Truncation of the Gaussian error distribution at +/-
          *truncation* * :math:`\sigma`
        seed:   Seed for the random generator. It *seed* is :py:obj:`None`, the
          :py:obj:`numpy.random.RandomState` singleton is used. If seed is an
          :py:class:`int`, a new :py:class:`~numpy.random.RandomState` instance
          is used. If seed is already a :py:class:`~numpy.random.Generator`
          like :py:obj:`at.random.common` then that instance is used.

    Other keyword arguments specify the kind of errors to be assigned. Error
    definitions may be given in 2 forms:

    * **a tuple of values**: (*systematic*, *random*). The *systematic* value
      is applied unchanged to all the selected elements, the *random* value is
      taken as the standard deviation of the generated error distribution,
    * **a single value**: it acts as the *random* part of the previous form,
      and there is no systematic error.

    For :py:class:`.Monitor` errors:

    The beam positions returned at :py:class:`.Monitor` locations are computed
    by applying in order a tilt :math:`\theta`, gains :math:`g_x`, :math:`g_y`
    and offsets :math:`o_x`, :math:`o_y` resulting in:

    :math:`\begin{pmatrix} x_{monitor} \\ y_{monitor} \end{pmatrix} =
    \begin{pmatrix} o_x \\ o_y \end{pmatrix} +
    \begin{pmatrix} g_x & 0 \\ -0 & g_y \end{pmatrix}
    \begin{pmatrix} cos\theta & sin\theta \\
    -sin\theta & cos\theta \end{pmatrix}
    \begin{pmatrix} x_{beam} \\ y_{beam} \end{pmatrix}`

    Note that this corresponds to the effect of:

    #. rotating the monitor device by :math:`-\theta`,
    #. then applying the gains,
    #. and finally shifting the monitor device by :math:`-o_x` horizontally and
       :math:`-o_y` vertically

    Keyword Args:
        BPMOffset:
        BPMGain:
        BPMTilt:

    """
    elements = ring.get_elements(key)
    for attr in _BPM_ATTRS + _ERR_ATTRS:
        val = kwargs.pop(attr, None)
        if val is not None:
            if isinstance(val, tuple):
                syst, rand = val
                rand = np.atleast_2d(rand)
                syst = np.atleast_2d(syst)
            else:
                rand = np.atleast_2d(val)
                syst = np.zeros(rand.shape)
            rv = _truncated_randn(size=(len(elements), rand.shape[-1]),
                                  truncation=truncation, random_state=seed)
            try:
                vals = syst + rv*rand
            except Exception as exc:
                exc.args = ('Attribute {0} does not accept value {1}: {2}'.format(
                    attr, val, exc),)
                raise
            for el, val in zip(elements, vals):
                setattr(el, attr, val)


def _rotmat(theta):
    cs = np.cos(theta)
    sn = np.sin(theta)
    return np.array([[cs, -sn], [sn, cs]])


def _apply_bpm_orbit_error(ring, refpts, orbit):
    for e, o6 in zip(refpts_iterator(ring, refpts), orbit):
        o6 = o6.reshape((-1, 6))
        if hasattr(e, 'BPMTilt'):
            o6[:, [0, 2]] = _rotmat(e.BPMTilt) @ o6[:, [0, 2]]
        if hasattr(e, 'BPMGain'):
            o6[: ,[0, 2]] *= e.BPMGain
        if hasattr(e, 'BPMOffset'):
            o6[: ,[0, 2]] += e.BPMOffset


def _apply_bpm_track_error(ring, refpts, trajectory):
    for traj in trajectory.T:
        _apply_bpm_orbit_error(ring, refpts, traj)


def _apply_alignment_errors(ring, **kwargs):
    refpts = [(hasattr(e, 'ShiftErr') or hasattr(e, 'RotationErr'))
              for e in ring]
    ring = ring.replace(refpts)
    alldef = kwargs.pop('all', True)
    for e in ring[refpts]:
        shift = getattr(e, 'ShiftErr')
        rots = getattr(e, 'RotationsErr')
        if shift is not None:
            shiftx = kwargs.pop('shiftx', alldef)*shift[0]
            shifty = kwargs.pop('shifty', alldef)*shift[1]
            shift_elem(e, shiftx, shifty, relative=True)
        if rots is not None:
            tilt = kwargs.pop('tilt', alldef)*rots[0]
            pitch = kwargs.pop('pitch', alldef) * rots[1]
            yaw = kwargs.pop('yaw', alldef) * rots[2]
            rotate_elem(e, tilt=tilt, pitch=pitch,
                        yaw=yaw, relative=True)
    return ring


def _apply_field_errors(ring, **kwargs):
    def sanitize(e):
        mo = max(len(e.PolynomA), len(e.PolynomB))
        e.PolynomA = np.pad(e.PolynomA, mo - len(e.PolynomA))
        e.PolynomB = np.pad(e.PolynomB, mo - len(e.PolynomB))
        e.MaxOrder = mo - 1

    def get_pol(e, pname, index):
        perr = np.copy(getattr(e, pname + 'Err'))
        if index is not None:
            perr[range(len(perr)) != index] = 0.0
        le = sorted((getattr(e, pname), perr), key=len)
        pn = np.copy(le[1])
        pn[:len(le[0])] += le[0]
        return pn

    def set_polerr(ring, pname, index):
        refpts = [hasattr(e, pname) and hasattr(e, pname+'Err') for e in ring]
        rr = ring.replace(refpts)
        for e in rr[refpts]:
            setattr(e, pname, get_pol(e, pname, index))
            sanitize(e)
        return rr

    alldef = kwargs.pop('all', True)
    for pol in ['A', 'B']:
        if kwargs.pop('Polynom'+pol, alldef):
            index = kwargs.pop('Index'+pol, None)
            ring = set_polerr(ring, 'Polynom'+pol, index)
    return ring
    

def enable_errors(ring: Lattice, **kwargs):
    if not ring.has_errors:
        ring = _apply_field_errors(ring, **kwargs)
        ring = _apply_alignment_errors(ring, **kwargs)
        ring._has_errors = True
    return ring


def has_errors(ring):
    return getattr(ring, '_has_errors', False)


def get_mean_std_err(ring: Lattice, key, attr, index=0):
    vals = [np.atleast_1d(getattr(e, attr, 0.0))[index]
            for e in ring.get_elements(key)]
    return np.mean(vals), np.std(vals)


def find_orbit_err(ring: Lattice, refpts: Refpts = None, **kwargs):
    ring = enable_errors(ring)
    orbit0, orbit = find_orbit(ring, refpts=refpts, **kwargs)
    _apply_bpm_orbit_error(ring, refpts, orbit)
    return orbit0, orbit


def get_optics_err(ring: Lattice, refpts: Refpts = None, **kwargs):
    ring = enable_errors(ring)
    ld0, bd, ld = get_optics(ring, refpts=refpts, **kwargs)
    _apply_bpm_orbit_error(ring, refpts, ld.closed_orbit)
    return ld0, bd, ld


def lattice_pass_err(ring: Lattice, r_in, nturns: int = 1,
                     refpts: Refpts = None, **kwargs):
    ring = enable_errors(ring)
    rout = lattice_pass(ring, r_in, nturns, refpts, **kwargs)
    _apply_bpm_track_error(ring, refpts, rout)
    return rout


Lattice.find_orbit_err = find_orbit_err
Lattice.get_optics_err = get_optics_err
Lattice.enable_errors = enable_errors
Lattice.has_errors = property(has_errors)
