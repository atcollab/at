# noinspection PyUnresolvedReferences
r"""Functions assigning errors to the magnets and monitors of a lattice,
enabling the errors and computing linear optics of the lattice with errors.

The processing of errors in PyAT is done in two steps:

#. The definition of errors: errors are stored in dedicated attributes of the
   lattice elements by :py:func:`assign_errors`.
   Nothing changes in tracking at this stage because the error attributes are
   still ignored,
#. The activation of errors by :py:func:`enable_errors`: Random errors are
   generated and are cumulated with systematic errors. Errors are stored in
   the standard element attributes of a new :py:class:`.Lattice`.

All error definitions may be given to the :py:func:`assign_errors` in 2
possible forms:

* **a tuple of values**: (*systematic*, *random*). The *systematic* value
  is applied unchanged to all the selected elements, the *random* value is
  taken as the standard deviation of the generated error distribution,
* **a single value**: it acts as the *random* part of the previous form,
  and there is no systematic error.

.. _monitor-errors:

Monitor errors
==============

The following errors are applied in that order:

BPMOffset:
    *[offset_x, offset_y]* array [m], all errors, electrical or mechanical,
    acting as if the monitor was shifted from the reference axis.
BPMTilt:
    *tilt* [rd], all errors, electrical or mechanical, acting as if the monitor
    was rotated around the z axis.
BPMGain:
    *[gain_x, gain_y]* array giving the magnification introduced by the monitor

This results in the following computation:

.. math::

    \begin{pmatrix} x_{monitor} \\ y_{monitor} \end{pmatrix} =
    \begin{pmatrix} g_x & 0 \\ 0 & g_y \end{pmatrix}
    \begin{pmatrix} cos\theta & sin\theta \\
    -sin\theta & cos\theta \end{pmatrix}
    \left[
    \begin{pmatrix} x_{beam} \\ y_{beam} \end{pmatrix} -
    \begin{pmatrix} o_x \\ o_y \end{pmatrix}
    \right]

.. note::
    Monitor errors may be applied to any element: when enabling errors, the
    beam position reported at the element location, whatever it is, will be
    affected by these attributes.

.. _magnet-alignment-errors:

Magnet alignment errors
=======================

ShiftErr:
    *[shift_x, shift_y]* array [m], magnet shift from the
    reference axis
RotationErr:
    *tilt* or *[tilt, pitch, yaw]* rotation around the z, x, y axis

.. _magnet-field-errors:

Magnet field errors
===================

PolynomBErr:
    Absolute field error. *PolynomAErr* is added to the magnet's *PolynomB*
PolynomAErr:
    Absolute field error. *PolynomBErr* is added to the magnet's *PolynomA*
ScalingPolynomBErr:
    Relative field error. *ScalingPolynomBErr* is the multipolar field error
    relative to the main magnet field component. When
    :py:func:`enabling the errors<enable_errors>`, *ScalingPolynomBErr* is
    multiplied by the magnet strengh and added to *PolynomB*
ScalingPolynomAErr:
    Relative field error. *ScalingPolynomAErr* is the skew multipolar field
    error relative to the main magnet field component.

.. admonition:: Example of use

    Let's imagine a standard quadrupole and see how error attributes can be
    related to the results of magnetic measurements:

    Systematic errors:
        The magnetic measurements at the nominal current give systematic
        octupole and 12-pole components in addition to the quadrupole.
        *PolynomB* looks like *[0, B2, 0, B4, 0, B6]*\ . The 1st step is to
        normalise with respect to the quadrupole component:

        *[0, 1, 0, B4/B2, 0, B6/B2]*

        The 1 goes to the main field, the rest goes to errors:

        >>> qp.PolynomB = [0, B2]
        >>> qp.ScalingPolynomBErr = (0, [0, 0, 0, B4/B2, 0, B6/B2])

        Note that the systematic error component must be in the 2nd position in
        the error tuple.

        Now assume that after tuning we set the quadrupole to C2 instead B2:

        >>> qp.PolynomB = [0, C2]

        :py:func:`enable_errors` will compute:

        *[0, C2] + C2*[0, 0, 0, B4/B2, 0, B6/B2] =
        [0, C2, 0, C2/B2*B4, 0, C2/B2*B6]*

        which is the exact scaling of the measured field by a ratio C2/B2.

    Random errors:
        In addition, we know that we have a random relative field integral error
        of 10\ :sup:`-3` and a random sextupole B3.  All this also scales with
        the magnet current, so we add a random component:

        >>> qp.ScalingPolynomBErr = ([0, 1.E-3, B3/B2, 1.E-3/B2, 0, 1.E-3/B2],
        ...                          [0, 0, 0, B4/B2, 0, B6/B2])

        After generating random numbers and multiplying by B2, we get the
        correct random contribution.

        Remark: in this example, if the field integral error results from
        mechanical magnet length, the random contributions to B2,
        B4 and B6 should be correlated, while the one of B3 is uncorrelated.
        Here a random value is generated for each component, so there is no
        correlation.

    Case of an arbitrary multipole magnet:
        Now let's consider a combined-function magnet, let's say a
        quadrupole-sextupole:

        >>> qp.PolynomB = [0, B2, B3]

        For errors, we can arbitrarily decide it's a :py:class:`.Quadrupole`
        and normalise the errors with B2:

        >>> qp.ScalingPolynomBErr = (0, [0, 0, 0, B4/B2, 0, B6/B2])

        which will then be multiplied by C2, strength of the quadrupole.

        or decide it's a :py:class:`.Sextupole` and normalise the errors with B3

        >>> qp.ScalingPolynomBErr = (0, [0, 0, 0, B4/B3, 0, B6/B3])

        which will be multiplied by C3, strength of the sextupole.

        Both choices are equivalent: the result will be **exactly the same**
        (provided the user correctly scales the main field *[0, C2, C3]*\ ).
        But we had to decide what is the "main order" of the magnet.


    Here we do not use the static error attribute. *ScalingPolynomBErr* is meant
    to be used for "single knob" magnets: the field scales with a single
    parameter, the magnet current. *ScalingPolynomBErr* ensures the correct
    scaling of errors without intervention, but it requires the "knob" to be
    associated with a scalar value: the strength of the main component,
    determined by the element class. For the :py:class:`.Multipole` class,
    there is no default order defined. The user has to set the
    *DefaultOrder* attribute to the order of its choice (dipole=0,
    quadrupole=1,…) and normalise the errors with this.
"""
import numpy as np
from typing import Optional, Union
from at.lattice import refpts_iterator, Lattice, Element, Refpts, elem_generator
from at.lattice import shift_elem, rotate_elem
from scipy.stats import truncnorm, norm


__all__ = ['apply_bpm_orbit_errors', 'apply_bpm_track_errors',
           'assign_errors', 'enable_errors']

_BPM_ATTRS = {'BPMGain': (2,), 'BPMOffset': (2,), 'BPMTilt': (1,)}
_ERR_ATTRS = {'PolynomBErr': None, 'PolynomAErr': None,
              'ScalingPolynomBErr': None, 'ScalingPolynomAErr': None,
              'ShiftErr': (2,), 'RotationErr': None}
_ALL_ATTRS = dict(**_BPM_ATTRS, **_ERR_ATTRS)

_SEL_ARGS = ('all', 'shiftx', 'shifty', 'tilt', 'pitch', 'yaw',
             'PolynomB', 'PolynomA', 'IndexB', 'IndexA',
             'BPMOffset', 'BPMGain', 'BPMTilt')


def _truncated_randn(truncation=None, **kwargs):
    if truncation is not None:
        return truncnorm.rvs(-truncation, truncation, **kwargs)
    else:
        return norm.rvs(**kwargs)


def assign_errors(ring: Lattice, refpts: Refpts, **kwargs):
    # noinspection PyUnresolvedReferences
    r"""Assign errors to selected elements

    The errors are stored as additional attributes but are not enabled for
    tracking (see :py:func:`enable_errors`).

    Args:
        ring:       Lattice description.
        refpts:     Observation points.
          See ":ref:`Selecting elements in a lattice <refpts>`"

    Other keyword arguments specify the kind of errors to be assigned. Error
    definitions may be given in 2 forms:

    * **a tuple of values**: (*systematic*, *random*). The *systematic* value
      is applied unchanged to all the selected elements, the *random* value is
      taken as the standard deviation of the generated error distribution,
    * **a single value**: it acts as the *random* part of the previous form,
      and there is no systematic error.

    .. rubric:: Monitor errors:

    See :ref:`monitor-errors` for a definition of monitor errors.

    Keyword Args:
        BPMOffset:  Offsets may be:

          * :py:class:`float` or (1,) float array: the same offset is applied
            in both planes to all elements
          * (2,) float array: [offset_x, offset_y] applied to all elements
          * (nelems, 1) or (nelems, 2) array: assign one value per element
            (this may be useful for systematic errors)
        BPMGain:    Monitor gain (scaling factor)

          * :py:class:`float` or (1,) float array: the same gain is applied
            in both planes to all elements
          * (2,) float array: [gain_x, gain_y] applied to all elements
          * (nelems, 1) or (nelems, 2) array: assign one value per element
            (this may be useful for systematic errors)
        BPMTilt:    :py:class:`float`, (1,) or (nelems, 1) float array

    .. rubric:: Magnet errors:

    See :ref:`magnet-alignment-errors` and :ref:`magnet-field-errors` for
    a definition of monitor errors.

    Keyword Args:
        ShiftErr:       Magnet shift [m]. Maybe:

          * :py:class:`float` or (1,) float array: the same shift is applied
            in both planes to all elements
          * (2,) float array: [shift\ :sub:`x`, shift\ :sub:`y`\ ] applied to
            all elements
          * (nelems, 1) or (nelems, 2) array: assign one value per element
            (this may be useful for systematic errors)
        RotationErr:    Magnet rotation

          * :py:class:`float` or (1,) float array: tilt error, no pitch nor yaw
          * (3,) float array: [tilt, pitch, yaw] rotations applied to all
            elements
          * (nelems, 1) or (nelems, 3) array: assign one value per element
            (this may be useful for systematic errors)
        PolynomBErr:    Absolute field error. *PolynomAErr* is added to the
          magnet's *PolynomB*
        PolynomAErr:    Absolute field error. *PolynomBErr* is added to the
          magnet's *PolynomA*
        ScalingPolynomBErr:    Field error. *ScalingPolynomBErr* is the
          multipolar field error relative to the main magnet field component.
          When :py:func:`enabling the errors<enable_errors>`,
          *ScalingPolynomBErr* is multiplied by the magnet strengh and added tp
          *PolynomB*
        ScalingPolynomAErr:    Field error. *ScalingPolynomAErr* is the skew
          multipolar field error relative to the main magnet field component.

    Examples:

        >>> qpoles =ring.get_cells(checktype(Quadrupole))
        >>> shift_err = 0.0001      # 100 um random position error (both planes)
        >>> rotation_err = 0.0002   # 200 urad random tilt error (both planes)
        >>> syst_PBErr = [0, 0, 0, 98]  # systematic octupole component
        >>> rand_PBErr = [0, 0.001]     # 1.e-3 relative gradient error
        >>> assign_errors(ring, qpoles,
        ... ShiftErr=shift_err,
        ... RotationErr=rotation_err,
        ... ScalingPolynomBErr=(syst_PBErr, rand_PBErr))



    See also:
        :py:func:`enable_errors`, :py:func:`get_optics_err`
    """
    elements = ring[refpts]
    nelems = len(elements)
    for attr, sz in _ALL_ATTRS.items():
        val = kwargs.pop(attr, None)
        if val is not None:
            if isinstance(val, tuple):
                syst, rand = val
                rand = np.atleast_2d(rand)
                syst = np.atleast_2d(syst)
            else:
                rand = np.atleast_2d(val)
                syst = np.zeros(rand.shape)
            if sz is None:
                szsyst = syst.shape[1:]
                szrand = rand.shape[1:]
            else:
                szsyst = szrand = sz
            syst = np.broadcast_to(syst, ((nelems,) + szsyst))
            rand = np.broadcast_to(rand, ((nelems,) + szrand))
            for el, s, r in zip(elements, syst, rand):
                setattr(el, attr, (s, r))


# noinspection PyProtectedMember
def apply_bpm_orbit_errors(ring: Lattice, refpts: Refpts, orbit):
    """Apply the Monitor errors on an orbit array

    Args:
        ring:           Lattice description.
        refpts:         Observation points.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        orbit:          (Nrefs, 6) closed orbit vector at each location
                        specified in *refpts*
    """

    def _rotmat(theta):
        cs = np.cos(theta)
        sn = np.sin(theta)
        return np.array([[cs, sn], [-sn, cs]])

    for elem, o6 in zip(refpts_iterator(ring, refpts), orbit):
        o6 = o6.reshape((-1, 6))
        if hasattr(elem, '_BPMOffset'):
            o6[:, [0, 2]] -= elem._BPMOffset
        if hasattr(elem, '_BPMTilt'):
            o6[:, [0, 2]] = o6[:, [0, 2]] @ _rotmat(elem._BPMTilt).T
        if hasattr(elem, '_BPMGain'):
            o6[:, [0, 2]] *= elem._BPMGain


def apply_bpm_track_errors(ring: Lattice, refpts: Refpts, trajectory):
    r"""Apply the Monitor errors on a trajectory array

    Args:
        ring:           Lattice description.
        refpts:         Observation points.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        trajectory:     (6, N, R, T) array containing output coordinates of
          N particles at R reference points for T turns.
    """
    # Loop on turns
    for traj in trajectory.T:
        apply_bpm_orbit_errors(ring, refpts, traj)


def _sysrand(value,
             truncation: float,
             seed: Union[int, np.random.Generator]):
    syst, rand = value
    rv = _truncated_randn(size=rand.shape,
                          truncation=truncation, random_state=seed)
    return syst + rv*rand


def _set_bpm_errors(elem: Element, errors,
                    truncation: float,
                    seed: Union[int, np.random.Generator],
                    **kwargs) -> Element:
    """Apply Monitor errors"""
    def set_locattr(el, err, inattr, outattr):
        if err is not None:
            delattr(el, inattr)
            value = _sysrand(err, truncation, seed)
            if kwargs.pop(inattr, alldef):
                setattr(el, outattr, value)

    alldef = kwargs.pop('all', True)
    offseterr, gainerr, tilterr = errors
    set_locattr(elem, offseterr, 'BPMOffset', '_BPMOffset')
    set_locattr(elem, gainerr, 'BPMGain', '_BPMGain')
    set_locattr(elem, tilterr, 'BPMTilt', '_BPMTilt')
    return elem


def _set_alignment_errors(elem: Element, errors,
                          truncation: float,
                          seed: Union[int, np.random.Generator],
                          **kwargs) -> Element:
    """Apply the defined magnet alignmemt errors"""
    alldef = kwargs.pop('all', True)
    shift, rots = errors

    if shift is not None:
        delattr(elem, 'ShiftErr')
        shift = _sysrand(shift, truncation, seed)
        shiftx = shift[0] + kwargs.pop('shiftx', alldef)
        shifty = shift[1] * kwargs.pop('shifty', alldef)
        shift_elem(elem, shiftx, shifty, relative=True)

    if rots is not None:
        delattr(elem, 'RotationErr')
        rots = _sysrand(rots, truncation, seed)
        tilt = rots[0] * kwargs.pop('tilt', alldef)
        pitch = rots[1] * kwargs.pop('pitch', alldef)
        yaw = rots[2] * kwargs.pop('yaw', alldef)
        rotate_elem(elem, tilt=tilt, pitch=pitch, yaw=yaw, relative=True)
    return elem


def _set_field_errors(elem: Element, errors_a, errors_b,
                      truncation: float,
                      seed: Union[int, np.random.Generator],
                      **kwargs):
    """Apply the defined field errors"""

    def sysrand(el: Element, error, attrname: str, enabled: bool, scale: float,
                index: Union[int, None], plist: list):
        """Append the errors attributes to the polynomial list"""
        def vmask(v, idx, pl):
            if idx is not None:
                if idx < len(v):
                    t = np.zeros(v.shape)
                    t[idx] = v[idx]
                    pl.append(t)
            elif len(v) > 0:
                pl.append(v)

        if error is not None:
            delattr(el, attrname)
            syst, rand = error
            rand = rand * _truncated_randn(size=rand.shape,
                                           truncation=truncation,
                                           random_state=seed)
            if enabled:
                vmask(scale*syst, index, plist)
                vmask(scale*rand, index, plist)

    def get_pol(el: Element, errs, pname: str, iname: str, **kwargs) -> list:
        """Build a list of all 5 polynomials for A or B"""
        errstatic, errdynamic = errs
        pstatic = pname+'Err'
        pdynamic = 'Scaling' + pstatic
        alldef = kwargs.pop('all', True)
        index = kwargs.pop(iname, None)
        enabled = kwargs.pop(pname, alldef)
        plist = [getattr(el, pname)]
        sysrand(el, errstatic, pstatic, enabled, 1.0, index, plist)
        sysrand(el, errdynamic, pdynamic, enabled, el.strength, index, plist)
        return plist

    def set_pol(el: Element, pname: str, plist: list, sz: int) -> None:
        """Sum up all polynomials and set the PolynomA/B attribute"""
        pn = np.zeros(sz)
        for plnm in plist:
            pn[:len(plnm)] += plnm
        setattr(el, pname, pn)

    alist = get_pol(elem, errors_a, 'PolynomA', 'IndexA', **kwargs)
    blist = get_pol(elem, errors_b, 'PolynomB', 'IndexB', **kwargs)
    psize = max(max(len(p) for p in alist), max(len(p) for p in blist))
    set_pol(elem, 'PolynomA', alist, psize)
    set_pol(elem, 'PolynomB', blist, psize)
    return elem


def enable_errors(ring: Lattice,
                  truncation: Optional[float] = None,
                  seed: Optional[Union[int, np.random.Generator]] = None,
                  **kwargs):
    r"""Enable magnet errors

    Returns a shallow copy of *ring* where magnet errors defined by
    :py:func:`assign_errors` are introduced in the magnet attributes.

    By default, all assigned errors are enabled, but keyword arguments allow a
    selection of the enabled errors.

    Args:
        ring:   Lattice description.
        truncation:     Truncation of the normal error distribution at +/-
          *truncation* * :math:`\sigma`. If None, no truncation is done.
        seed:   Seed for the random generator. It *seed* is :py:obj:`None`, a
          new :py:class:`~numpy.random.Generator` instance with fresh,
           unpredictable entropy is created and used. If seed is an
          :py:class:`int`, a new :py:class:`~numpy.random.Generator` instance
          with its initial state given by *seed*. If seed is already a
          :py:class:`~numpy.random.Generator` like :py:obj:`at.random.common`
          then that instance is used.

    Keyword Args:
        all (bool):         Set the default value for all the specific error
          flags. Default: :py:obj:`True`
        shiftx (bool):      enable horizontal shift errors. Default: *all*
        shifty (bool):      enable vertical shift errors. Default: *all*
        tilt (bool):        enable tilt errors. Default: *all*
        pitch (bool):       enable pitch errors. Default: *all*
        yaw (bool):         enable yaw errors. Default: *all*
        PolynomB (bool):    enable polynomial errors. Default: *all*
        IndexB (Optional[int]): restrict the polynomial errors to the
          specified index. Default: :py:obj:`None` meaning the whole polynom
        PolynomA (bool):    enable polynomial errors. Default: *all*
        IndexA (Optional[int]): restrict the skew polynomial error to the
          specified index. Default: :py:obj:`None` meaning the whole polynom
        BPMOffset (bool):   enable Monitor offset errors. Default: *all*
        BPMGain (bool):     enable Monitor gain errors. Default: *all*
        BPMTilt (bool):     enable Monitor tilt errors. Default: *all*

    Examples:
        >>> ringerr = enable_errors(ring)

        All errors are enabled in *ringerr*

        >>> ringerr = enable_errors(all=False, ShiftErr=True)

        Only magnet shift errors are active

        >>> ringerr = enable_errors(all=False)

        No magnet errors are enabled. Only monitor errors will be active.

    See also:
        :py:func:`assign_errors`, :py:func:`get_optics_err`
    """
    def error_generator(trunc: float,
                        sd: Union[int, np.random.Generator]
                        ):
        def error_attributes(el, *args):
            return tuple(getattr(el, attr, None) for attr in args)

        for elem in ring:
            pola_errors = error_attributes(elem,
                                           'PolynomAErr', 'ScalingPolynomAErr')
            polb_errors = error_attributes(elem,
                                           'PolynomBErr', 'ScalingPolynomBErr')
            if not all(a is None for a in pola_errors + polb_errors):
                elem = _set_field_errors(elem.deepcopy(),
                                         pola_errors, polb_errors,
                                         truncation=trunc,
                                         seed=sd, **kwargs)
            align_errors = error_attributes(elem, 'ShiftErr', 'RotationErr')
            if not all(a is None for a in align_errors):
                elem = _set_alignment_errors(elem.deepcopy(), align_errors,
                                             truncation=trunc,
                                             seed=sd, **kwargs)
            bpm_errors = error_attributes(elem,
                                          'BPMOffset', 'BPMGain', 'BPMTilt')
            if not all(a is None for a in bpm_errors):
                elem = _set_bpm_errors(elem.deepcopy(), bpm_errors,
                                       truncation=trunc,
                                       seed=sd, **kwargs)
            yield elem

    if not ring.has_errors:
        gen = np.random.default_rng(seed)
        ring = Lattice(elem_generator,
                       error_generator(truncation, gen),
                       iterator=ring.attrs_filter, _has_errors=True)
    return ring


def has_errors(ring):
    """Get the error status of a lattice"""
    return getattr(ring, '_has_errors', False)


def get_mean_std_err(ring: Lattice, key, attr, index=0):
    vals = [np.atleast_1d(getattr(e, attr, 0.0))[index]
            for e in ring.select(key)]
    return np.mean(vals), np.std(vals)


def _sort_flags(kwargs):
    errargs = {}
    for key in _SEL_ARGS:
        if key in kwargs:
            errargs[key] = kwargs.pop(key)
    return kwargs, errargs


Lattice.assign_errors = assign_errors
Lattice.enable_errors = enable_errors
Lattice.has_errors = property(has_errors, doc="Error status")
