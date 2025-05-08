import numpy
from .lattice_object import Lattice
from .elements import Element
from .utils import Refpts, BoolRefpts, Uint32Refpts
from .utils import checkattr, get_uint32_index, get_bool_index
from typing import Sequence

__all__ = ['get_cells', 'get_refpts']


# noinspection PyIncorrectDocstring
def get_cells(ring: Sequence[Element], refpts: Refpts, *args,
              regex=False) -> BoolRefpts:
    # noinspection PyShadowingNames
    r"""
    get_cells(ring, filtfunc) -> BoolRefpts
    get_cells(ring, element_type) -> BoolRefpts
    get_cells(ring, attrname) -> BoolRefpts
    get_cells(ring, attrname, attrvalue) -> BoolRefpts
    Returns a bool array of element indices, selecting ring elements.

    Deprecated: :pycode:`get_cells(ring, refpts)` is
    :pycode:`ring.bool_refpts(refpts)` except for :py:obj:`str` arguments:
    :pycode:`get_cells(ring, attrname [, attrvalue])` is
    :pycode:`ring.bool_refpts(checkattr(strkey [, attrvalue]))`

    Parameters:
        ring (Sequence[Element]):       Lattice description
        filtfunc (ElementFilter):       Filter function. Selects
          :py:class:`.Element`\ s satisfying the filter function
        element_type (Type[Element]):   Element type
        attrname (str):                 Attribute name
        attrvalue (Any):                Attribute value. If absent, select the
          presence of an *attrname* attribute. If present, select
          :py:class:`.Element`\ s with :pycode:`attrname == attrvalue`.
        regex: Use regular expression for refpts string matching;
            Default: False (Unix shell-style wildcards)

    Returns:
        bool_refs (BoolRefpts):  numpy Array of :py:obj:`bool` with length
          len(ring)+1

    Examples:

        >>> refpts = get_cells(ring, 'Frequency')

        Returns a numpy array of booleans where all elements having a
        :pycode:`Frequency` attribute :py:obj:`True`

        >>> refpts = get_cells(ring, 'K', 0.0)

        Returns a numpy array of booleans where all elements having a
        :pycode:`K` attribute equal to 0.0 are :py:obj:`True`

    See also:
        :py:meth:`.Lattice.bool_refpts`, :py:meth:`.Lattice.uint32_refpts`
    """
    if isinstance(refpts, str):
        refpts = checkattr(refpts, *args)
    return get_bool_index(ring, refpts, regex=regex)


# noinspection PyUnusedLocal,PyIncorrectDocstring
def get_refpts(ring: Sequence[Element], refpts: Refpts,
               regex=False) -> Uint32Refpts:
    r"""Return a :py:obj:`~numpy.uint32` array of element indices selecting
    ring elements.

    Deprecated: :pycode:`get_elements(ring, refpts)` is
    :pycode:`ring.uint32_refpts(refpts)`

    Parameters:
        ring:           Lattice description
        refpts:         Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        regex: Use regular expression for refpts string matching;
            Default: False (Unix shell-style wildcards)

    Returns:
        uint32_refs (Uint32Refs):    :py:obj:`~numpy.uint32` numpy array as
          long as the number of refpts

    See also:
        :py:meth:`.Lattice.uint32_refpts`, :py:meth:`.Lattice.bool_refpts`
    """
    return get_uint32_index(ring, refpts, regex=regex)


def rotate_elem(elem: Element, tilt: float = 0.0, pitch: float = 0.0,
                yaw: float = 0.0, relative: bool = False) -> None:
    r"""Set the tilt, pitch, and yaw angle of an :py:class:`.Element`.
    The tilt is a rotation around the *s*-axis, the pitch is a
    rotation around the *x*-axis, and the yaw is a rotation around
    the *y*-axis.

    A positive angle represents a clockwise rotation when
    looking in the direction of the rotation axis.

    The transformations are not all commutative, the pitch and yaw
    are applied first, and the tilt is always the last transformation
    applied. The element is rotated around its mid-point.

    If *relative* is :py:obj:`True`, the previous angle and shifts
    are rebuilt from the *R* and *T* matrix and incremented by the
    input arguments.

    The shift is always conserved regardless of the value of *relative*.

    The transformations are applied by changing the particle coordinates
    at the entrance of the element and restoring them at the end. Following
    the small angles approximation, the longitudinal shift of the particle
    coordinates is neglected and the element length is unchanged.

    Parameters:
        elem:           Element to be tilted
        tilt:           Tilt angle [rad]
        pitch:          Pitch angle [rad]
        yaw:            Yaw angle [rad]
        relative:       If :py:obj:`True`, the rotation is added to the
          previous one
    """
    # noinspection PyShadowingNames
    def _get_rm_tv(le, tilt, pitch, yaw):
        tilt = numpy.around(tilt, decimals=15)
        pitch = numpy.around(pitch, decimals=15)
        yaw = numpy.around(yaw, decimals=15)
        ct, st = numpy.cos(tilt), numpy.sin(tilt)
        ap, ay = 0.5*le*numpy.tan(pitch), 0.5*le*numpy.tan(yaw)
        rr1 = numpy.asfortranarray(numpy.diag([ct, ct, ct, ct, 1.0, 1.0]))
        rr1[0, 2] = st
        rr1[1, 3] = st
        rr1[2, 0] = -st
        rr1[3, 1] = -st
        rr2 = rr1.T
        t1 = numpy.array([ay, numpy.sin(-yaw), -ap, numpy.sin(pitch), 0, 0])
        t2 = numpy.array([ay, numpy.sin(yaw), -ap, numpy.sin(-pitch), 0, 0])
        rt1 = numpy.eye(6, order='F')
        rt1[1, 4] = t1[1]
        rt1[3, 4] = t1[3]
        rt2 = numpy.eye(6, order='F')
        rt2[1, 4] = t2[1]
        rt2[3, 4] = t2[3]
        return rr1 @ rt1, rt2 @ rr2, t1, t2

    tilt0 = 0.0
    pitch0 = 0.0
    yaw0 = 0.0
    t10 = numpy.zeros(6)
    t20 = numpy.zeros(6)
    if hasattr(elem, 'R1') and hasattr(elem, 'R2'):
        rr10 = numpy.eye(6, order='F')
        rr10[:4, :4] = elem.R1[:4, :4]
        rt10 = rr10.T @ elem.R1
        tilt0 = numpy.arctan2(rr10[0, 2], rr10[0, 0])
        yaw0 = numpy.arcsin(-rt10[1, 4])
        pitch0 = numpy.arcsin(rt10[3, 4])
        _, _, t10, t20 = _get_rm_tv(elem.Length, tilt0, pitch0, yaw0)
    if hasattr(elem, 'T1') and hasattr(elem, 'T2'):
        t10 = elem.T1-t10
        t20 = elem.T2-t20
    if relative:
        tilt += tilt0
        pitch += pitch0
        yaw += yaw0

    r1, r2, t1, t2 = _get_rm_tv(elem.Length, tilt, pitch, yaw)
    elem.R1 = r1
    elem.R2 = r2
    elem.T1 = t1+t10
    elem.T2 = t2+t20


Lattice.uint32_refpts = get_uint32_index
Lattice.bool_refpts = get_bool_index
Lattice.get_cells = get_cells
Lattice.get_refpts = get_refpts
