"""Element translations and rotations.

.. Caution::

    The geometrical transformations are not commutative. When combining several
    transformations on the same element by using multiple function calls, the
    transformations are applied in a fixed order, independent of the order of the
    function calls:

    .. centered:: translations -> tilt (z-axis) -> yaw (y-axis) -> pitch (x-axis)
"""

from __future__ import annotations

__all__ = [
    "ReferencePoint",
    "transform_elem",
    "shift_elem",
    "tilt_elem",
    "set_shift",
    "set_tilt",
    "set_rotation",
    "transform_options",
]

from enum import Enum
from collections.abc import Sequence

import numpy as np

from .elements import Element
from .utils import Refpts, All, refpts_iterator, _refcount

_x_axis = np.array([1.0, 0.0, 0.0])
_y_axis = np.array([0.0, 1.0, 0.0])
_z_axis = np.array([0.0, 0.0, 1.0])


class ReferencePoint(Enum):
    """Definition of the reference point for the geometric transformations."""
    CENTRE = 0  #: Origin at the centre of the element.
    ENTRANCE = 1  #: Origin at the entrance of the element.


class _TransFormOptions:
    rounding = None
    referencepoint = ReferencePoint.CENTRE


transform_options = _TransFormOptions()


# noinspection PyPep8Naming
def _rotation(rotations):
    """
    The implementation follows the one described in:
    https://doi.org/10.1016/j.nima.2022.167487
    All the comments featuring 'Eq' points to the paper's equations.

    3D rotation matrix (extrinsic rotation) using the Tait–Bryan angles
    convention.
    For more details, refer to https://en.wikipedia.org/wiki/Euler_angles
    Corresponds to Eq. (3)
    alpha: Rotation about the X-axis (pitch).
    beta: Rotation about the Y-axis (yaw).
    gamma: Rotation about the Z-axis (roll/tilt).
    """
    alpha, beta, gamma = rotations  # ZYX intrinsic rotations (pitch, yaw, tilt)
    R_x = np.array(
        [
            [1, 0, 0],
            [0, np.cos(alpha), -np.sin(alpha)],
            [0, np.sin(alpha), np.cos(alpha)],
        ]
    )
    R_y = np.array(
        [[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]]
    )
    R_z = np.array(
        [
            [np.cos(gamma), -np.sin(gamma), 0],
            [np.sin(gamma), np.cos(gamma), 0],
            [0, 0, 1],
        ]
    )
    return R_x @ R_y @ R_z


# noinspection PyPep8Naming
def _translation_vector(ld, r3d, offsets, X_axis, Y_axis):
    """
    The implementation follows the one described in:
    https://doi.org/10.1016/j.nima.2022.167487
    All the comments featuring 'Eq' points to the paper's equations.

    Translation vector resulting from the joint effect of a longitudinal
    displacement (in the rotated frame), 3D offsets, and the 3D rotation matrix.
    Corresponds to Eqs. (8-11)
    ld: Longitudinal displacement [m]
    r3d: 3D rotation matrix
    offsets: 3D offsets [m]
    X_axis: X unit axis in rotated frame expressed in the xyz coordinate system
    Y_axis: Y unit axis in rotated frame expressed in the xyz coordinate system
    """
    tD0 = np.array([-offsets @ X_axis, 0, -offsets @ Y_axis, 0, 0, 0])
    T0 = np.array(
        [
            ld * r3d[2, 0] / r3d[2, 2],
            r3d[2, 0],
            ld * r3d[2, 1] / r3d[2, 2],
            r3d[2, 1],
            0,
            ld / r3d[2, 2],
        ]
    )
    return T0 + tD0


def _r_matrix(ld, r3d):
    """
    The implementation follows the one described in:
    https://doi.org/10.1016/j.nima.2022.167487
    All the comments featuring 'Eq' points to the paper's equations.

    Rotation matrix operator (R1, R2).
    Can take into account the effect of a longitudinal displacement.
    ld: Longitudinal displacement [m]
    r3d: 3D rotation matrix
    Corresponds to Eq. (9)
    """
    return np.array(
        [
            [
                r3d[1, 1] / r3d[2, 2],
                ld * r3d[1, 1] / r3d[2, 2] ** 2,
                -r3d[0, 1] / r3d[2, 2],
                -ld * r3d[0, 1] / r3d[2, 2] ** 2,
                0,
                0,
            ],
            [0, r3d[0, 0], 0, r3d[1, 0], r3d[2, 0], 0],
            [
                -r3d[1, 0] / r3d[2, 2],
                -ld * r3d[1, 0] / r3d[2, 2] ** 2,
                r3d[0, 0] / r3d[2, 2],
                ld * r3d[0, 0] / r3d[2, 2] ** 2,
                0,
                0,
            ],
            [0, r3d[0, 1], 0, r3d[1, 1], r3d[2, 1], 0],
            [0, 0, 0, 0, 1, 0],
            [
                -r3d[0, 2] / r3d[2, 2],
                -ld * r3d[0, 2] / r3d[2, 2] ** 2,
                -r3d[1, 2] / r3d[2, 2],
                -ld * r3d[1, 2] / r3d[2, 2] ** 2,
                0,
                1,
            ],
        ]
    )


def _tilt_frame_mat(rots: float) -> None:
        cs = np.cos(rots)
        sn = np.sin(rots)
        rm = np.asfortranarray(np.diag([cs, cs, cs, cs, 1.0, 1.0]))
        rm[0, 2] = sn
        rm[1, 3] = sn
        rm[2, 0] = -sn
        rm[3, 1] = -sn
        return rm
        

# noinspection PyPep8Naming
def transform_elem(
    elem: Element,
    dx: float | None = None,
    dy: float | None = None,
    dz: float | None = None,
    tilt: float | None = None,
    pitch: float | None = None,
    yaw: float | None = None,
    tilt_frame: float | None = None,
    *,
    reference: ReferencePoint | None = None,
    relative: bool = False,
) -> None:
    r"""Set the translations and rotations of an :py:class:`.Element`.

    The tilt is a rotation around the *s*-axis, the pitch is a
    rotation around the *x*-axis, and the yaw is a rotation around the *y*-axis.

    A positive angle represents a clockwise rotation when
    looking in the direction of the rotation axis.

    The transformations are not all commutative. The translations are applied before
    the rotations. The rotations are applied in the order *Z* -> *Y* -> *X*
    (tilt -> yaw -> pitch). The element is rotated around its mid-point. The mid-point
    can either be the element entrance or its centre (axis joining the entry and exit
    points of the element).

    If *relative* is :py:obj:`True`, the previous translations and angles
    are incremented by the input arguments.

    PyAT describes the ultra-relativistic beam dynamics in 6D phase space
    coordinates, which differ from 3D spatial angles in an expansion with
    respect to the energy to first order by a factor :math:`(1 + \delta)`, where
    :math:`\delta` is the relative momentum deviation. This introduces
    spurious dispersion (angle proportional to :math:`\delta`).

    The implementation follows the one described in:
    https://doi.org/10.1016/j.nima.2022.167487
    All the comments featuring 'Eq' point to the paper's equations.

    Parameters:
        elem:           Element to be transformed.
                        Default: :py:obj:`ReferencePoint.CENTRE`
        dx:             Horizontal shift [m]. Default: no change.
        dy:             Vertical shift [m]. Default: no change.
        dz:             Longitudinal shift [m]. Default: no change.
        tilt:           Tilt angle [rad]. Default: no change.
        pitch:          Pitch angle [rad]. Default: no change.
        yaw:            Yaw angle [rad]. Default: no change
        tilt_frame:     Tilt angle of the reference frame [rad]. Useful
                        to generate vertical bending magnets for example
                        Default no change
        reference:      Transformation reference, either
                        :py:obj:`ReferencePoint.CENTRE` or
                        :py:obj:`ReferencePoint.ENTRANCE`.
        relative:       If :py:obj:`True`, the input values are added to the
          previous ones.

    .. Attention::

        When combining several transformations by using multiple calls to
        :py:func:`transform_elem`, the *reference* argument must be identical for all.
        Setting the *reference* will therefore affect all transformations for this element.
    """
    if reference is not None:
        elem.ReferencePoint = reference

    if relative:
        def _set(ini, val):
            return ini if val is None else ini + val
    else:
        def _set(ini, val):
            return ini if val is None else val
        
    if transform_options.rounding is not None:
        if dx is not None: dx = np.round(dx, transform_options.rounding)
        if dy is not None: dy = np.round(dy, transform_options.rounding)
        if dz is not None: dz = np.round(dz, transform_options.rounding)
        if pitch is not None: pitch = np.round(pitch, transform_options.rounding)
        if yaw is not None: yaw  = np.round(yaw, transform_options.rounding)
        if tilt is not None: tilt = np.round(tilt, transform_options.rounding)

    elem_length = getattr(elem, "Length", 0)
    elem_bending_angle = getattr(elem, "BendingAngle", 0)

    # Rotated transfer matrix linked to a bending angle
    RB = _rotation([0, -elem_bending_angle, 0])  # Eq. (12)
    RB_half = _rotation([0, -elem_bending_angle / 2, 0])  # Eq. (12)

    # Get the current transformation
    offsets0 = [elem.dx, elem.dy, elem.dz]
    tilt0 = elem.tilt
    yaw0 = elem.yaw
    pitch0 = elem.pitch

    # Apply new offsets and rotations (XYZ intrinsic order)
    offsets = np.array([_set(v0, v) for v0, v in zip(offsets0, [dx, dy, dz])])
    rotations = [_set(pitch0, pitch), _set(yaw0, yaw), _set(tilt0, tilt)]

    setattr(elem, "_dx", offsets[0])
    setattr(elem, "_dy", offsets[1])
    setattr(elem, "_dz", offsets[2])
    setattr(elem, "_pitch", rotations[0])
    setattr(elem, "_yaw", rotations[1])
    setattr(elem, "_tilt", rotations[2])

    if elem.ReferencePoint is ReferencePoint.CENTRE:
        # Compute entrance rotation matrix in the rotated frame
        r3d_entrance = RB_half @ _rotation(rotations) @ RB_half.T  # Eq. (31)

        if elem_bending_angle:
            Rc = elem_length / elem_bending_angle
            OO0 = Rc * np.sin(elem_bending_angle / 2) * RB_half @ _z_axis  # Eq. (34)
            P0P = (
                -Rc * np.sin(elem_bending_angle / 2) * r3d_entrance @ RB_half @ _z_axis
            )  # Eq. (36)
        else:
            OO0 = elem_length / 2 * _z_axis  # Eq. (34)
            P0P = -elem_length / 2 * r3d_entrance @ _z_axis  # Eq. (36)

        # Transform offset to magnet entrance
        OP = OO0 + P0P + RB_half @ offsets  # Eq. (33)

    elif elem.ReferencePoint is ReferencePoint.ENTRANCE:
        r3d_entrance = _rotation(rotations)  # Eq. (3)
        OP = offsets  # Eq. (2)
    else:
        raise ValueError(
            "Unsupported reference, please choose either "
            "ReferencePoint.CENTRE or "
            "ReferencePoint.ENTRANCE."
        )

    # R1, T1
    # XYZ - axes unit - vectors expressed in the xyz coordinate system
    X_axis = r3d_entrance @ _x_axis
    Y_axis = r3d_entrance @ _y_axis
    Z_axis = r3d_entrance @ _z_axis

    ld_entrance = Z_axis @ OP  # Eq. (33)

    R1 = _r_matrix(ld_entrance, r3d_entrance)
    T1 = np.linalg.inv(R1) @ _translation_vector(
        ld_entrance, r3d_entrance, OP, X_axis, Y_axis
    )

    # R2, T2
    # XYZ - axes unit - vectors expressed in the xyz coordinate system
    X_axis = RB @ _x_axis
    Y_axis = RB @ _y_axis
    Z_axis = RB @ _z_axis

    r3d_exit = RB.T @ r3d_entrance.T @ RB  # Eq. (18) or (32)

    if elem_bending_angle:
        Rc = elem_length / elem_bending_angle
        OPp = np.array(
            [
                Rc * (np.cos(elem_bending_angle) - 1),
                0,
                elem_length * np.sin(elem_bending_angle) / elem_bending_angle,
            ]
        )  # Eq. (24)
    else:
        OPp = elem_length * _z_axis  # Eq. (24)

    OOp = r3d_entrance @ OPp + OP  # Eq. (25)
    OpPp = OPp - OOp

    ld_exit = Z_axis @ OpPp  # Eq. (23) or (37)

    R2 = _r_matrix(ld_exit, r3d_exit)
    T2 = _translation_vector(ld_exit, r3d_exit, OpPp, X_axis, Y_axis)

    # Update element
    elem.R1 = R1
    elem.R2 = R2
    elem.T1 = T1
    elem.T2 = T2

    if tilt_frame is not None:
        tf_mat = _tilt_frame_mat(tilt_frame)
        elem.R1 = tf_mat @ elem.R1
        elem.R2 = elem.R2 @ tf_mat.T 


def tilt_elem(elem: Element, rots: float | None = None, rots_frame: float | None = None,
              relative: bool = False, reference: ReferencePoint | None = None) -> None:
    r"""Set the tilt angle :math:`\theta` of an :py:class:`.Element`

    The rotation matrices are stored in the :pycode:`R1` and :pycode:`R2`
    attributes.

    :math:`R_1=\begin{pmatrix} cos\theta & sin\theta \\
    -sin\theta & cos\theta \end{pmatrix}`,
    :math:`R_2=\begin{pmatrix} cos\theta & -sin\theta \\
    sin\theta & cos\theta \end{pmatrix}`

    Parameters:
        elem:           Element to be tilted
        rots:           Tilt angle :math:`\theta` [rd]. *rots* > 0 corresponds
          to a corkscrew rotation of the element looking in the direction of
          the beam. Use :py:obj:`None` to keep the current value.
        rots_frame:     Tilt angle :math:`\theta` [rd]. *rots* > 0 corresponds
          to a corkscrew rotation of the reference frame looking in the direction of
          the beam. Use :py:obj:`None` to keep the current value.
        relative:       If :py:obj:`True`, the rotation is added to the
          existing one
        reference:      Transformation reference, either
          :py:obj:`ReferencePoint.CENTRE` or
          :py:obj:`ReferencePoint.ENTRANCE`.

    See Also:
        :py:func:`shift_elem`
        :py:func:`.transform_elem`
    """
    transform_elem(elem, tilt=rots, tilt_frame=rots_frame, relative=relative, reference=reference)


def shift_elem(
    elem: Element,
    dx: float | None = None,
    dy: float | None = None,
    dz: float | None = None,
    *,
    reference: ReferencePoint | None = None,
    relative: bool = False,
) -> None:
    r"""Sets the translations of an :py:class:`.Element`

    The translation vectors are stored in the :pycode:`T1` and :pycode:`T2`
    attributes.

    Parameters:
        elem:           Element to be shifted
        dx:             Horizontal translation [m]. Default no change.
        dy:             Vertical translation [m]. Default no change.
        dz:             Longitudinal translation [m]. Default no change.
        reference:      Transformation reference, either
                        :py:obj:`ReferencePoint.CENTRE` or
                        :py:obj:`ReferencePoint.ENTRANCE`.
        relative:       If :py:obj:`True`, the translation is added to the
          existing one

    See Also:
        :py:func:`tilt_elem`
        :py:func:`.transform_elem`
    """
    transform_elem(elem, reference=reference, dx=dx, dy=dy, dz=dz, relative=relative)


def set_rotation(
    ring: Sequence[Element],
    tilts=None,
    pitches=None,
    yaws=None,
    tilts_frame=None,
    *,
    refpts: Refpts = All,
    reference: ReferencePoint | None = None,
    relative=False,
) -> None:
    r"""Sets the rotations of a list of elements.

    Parameters:
        ring:        Lattice description.
        tilts:       Scalar or Sequence of tilt values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        tilts_frame: Scalar or Sequence of reference frame tilt values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        pitches:     Scalar or Sequence of pitch values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        yaws:        Scalar or Sequence of yaw values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        refpts:      Element selection key.
          See ":ref: `Selecting elements in a lattice <refpts>`"
        reference:   Transformation reference, either
                     :py:obj:`ReferencePoint.CENTRE` or
                     :py:obj:`ReferencePoint.ENTRANCE`.
        relative:    If :py:obj:`True`, the rotations are added to the existing ones.

    .. versionadded:: 0.7.0
       The *refpts* argument

    See Also:
        :py:func:`set_tilt`
        :py:func:`set_shift`
    """
    nb = _refcount(ring, refpts, endpoint=False)
    tilts = np.broadcast_to(tilts, (nb,))
    pchs = np.broadcast_to(pitches, (nb,))
    yaws = np.broadcast_to(yaws, (nb,))
    tilts_frame = np.broadcast_to(tilts_frame, (nb,))
    for el, tilt, pitch, yaw, tilt_frame in zip(refpts_iterator(ring, refpts), tilts, pchs, yaws, tilts_frame):
        transform_elem(el, reference=reference, tilt=tilt, pitch=pitch, yaw=yaw,
                       tilt_frame=tilt_frame, relative=relative)


def set_tilt(
    ring: Sequence[Element], tilts: float | None = None, tilts_frame: float | None = None, *,
    refpts: Refpts = All, reference: ReferencePoint | None = None,
    relative=False
) -> None:
    r"""Sets the tilts of a list of elements.

    Parameters:
        ring:        Lattice description.
        tilts:       Scalar or Sequence of tilt values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        tilts_frame: Scalar or Sequence of reference frame tilt values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        refpts:      Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        reference:   Transformation reference, either
          :py:obj:`ReferencePoint.CENTRE` or
          :py:obj:`ReferencePoint.ENTRANCE`.
        relative:    If :py:obj:`True`, the rotation is added to the existing one.

    .. versionadded:: 0.7.0
       The *refpts* argument

    See Also:
        :py:func:`set_rotation`
        :py:func:`set_shift`
    """
    nb = _refcount(ring, refpts, endpoint=False)
    tilts = np.broadcast_to(tilts, (nb,))
    tilts_frame = np.broadcast_to(tilts_frame, (nb,))
    for el, tilt, tilt_frame in zip(refpts_iterator(ring, refpts), tilts, tilts_frame):
        transform_elem(el, reference=reference, tilt=tilt, tilt_frame=tilt_frame, relative=relative)


def set_shift(
    ring: Sequence[Element], dxs, dys, dzs=None, *,
    refpts: Refpts = All, reference: ReferencePoint | None = None, relative=False
) -> None:
    r"""Sets the translations of a list of elements.

    Parameters:
        ring:       Lattice description.
        dxs:        Scalar or Sequence of horizontal translations values applied
          to the selected elements. Use :py:obj:`None` to keep the current values [m].
        dys:        Scalar or Sequence of vertical translations values applied
          to the selected elements. Use :py:obj:`None` to keep the current values [m].
        dzs:        Scalar or Sequence of longitudinal translations values applied
          to the selected elements. Use :py:obj:`None` to keep the current values [m].
        refpts:     Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        reference:  Transformation reference, either
                    :py:obj:`ReferencePoint.CENTRE` or
                    :py:obj:`ReferencePoint.ENTRANCE`.
        relative:   If :py:obj:`True`, the translation is added to the
          existing one.

    .. versionadded:: 0.7.0
       The *refpts* argument

    See Also:
        :py:func:`set_rotation`
        :py:func:`set_tilt`
    """
    nb = _refcount(ring, refpts, endpoint=False)
    dxs = np.broadcast_to(dxs, (nb,))
    dys = np.broadcast_to(dys, (nb,))
    dzs = np.broadcast_to(dzs, (nb,))
    for el, dx, dy, dz in zip(refpts_iterator(ring, refpts), dxs, dys, dzs):
        transform_elem(el, reference=reference, dx=dx, dy=dy, dz=dz, relative=relative)


def _get_referencePoint(elem: Element) -> ReferencePoint: 
    "Rotation reference point"
    idx = getattr(elem, "_referencepoint", transform_options.referencepoint.value)
    return list(ReferencePoint)[idx]


def _set_referencePoint(elem: Element, value: ReferencePoint) -> None: 
    setattr(elem, "_referencepoint", value.value)


def _get_dx(elem: Element) -> float:
    """Horizontal element shift"""
    return getattr(elem, "_dx", 0.0)


def _set_dx(elem: Element, value: float) -> None:
    transform_elem(elem, dx=value)


def _get_dy(elem: Element) -> float:
    """Vertical element shift"""
    return getattr(elem, "_dy", 0.0)


def _set_dy(elem: Element, value: float) -> None:
    transform_elem(elem, dy=value)


def _get_dz(elem: Element) -> float:
    """Longitudinal element shift"""
    return getattr(elem, "_dz", 0.0)


def _set_dz(elem: Element, value: float) -> None:
    transform_elem(elem, dz=value)


def _get_tilt(elem: Element) -> float:
    """Element tilt"""
    return getattr(elem, "_tilt", 0.0)


def _set_tilt(elem: Element, value: float) -> None:
    transform_elem(elem, tilt=value)


def _get_pitch(elem: Element) -> float:
    """Element pitch"""
    return getattr(elem, "_pitch", 0.0)


def _set_pitch(elem: Element, value: float) -> None:
    transform_elem(elem, pitch=value)


def _get_yaw(elem: Element) -> float:
    """Element yaw"""
    return getattr(elem, "_yaw", 0.0)


def _set_yaw(elem: Element, value: float) -> None:
    transform_elem(elem, yaw=value)


def _get_tilt_frame(elem: Element) -> float:
    """Element tilt frame, different from tilt only for bends"""
    return getattr(elem, "_tilt_frame", 0.0)


def _set_tilt_frame(elem: Element, value: float) -> None:
    transform_elem(elem, tiltframe=value)


Element.transform = transform_elem
Element.ReferencePoint = property(_get_referencePoint, _set_referencePoint)
Element.dx = property(_get_dx, _set_dx)
Element.dy = property(_get_dy, _set_dy)
Element.dz = property(_get_dz, _set_dz)
Element.tilt = property(_get_tilt, _set_tilt)
Element.pitch = property(_get_pitch, _set_pitch)
Element.yaw = property(_get_yaw, _set_yaw)
Element.tilt_frame = property(_get_tilt_frame, _set_tilt_frame)