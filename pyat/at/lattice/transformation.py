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
    "get_offsets_rotations",
    "transform_elem",
    "shift_elem",
    "tilt_elem",
    "set_shift",
    "set_tilt",
    "set_rotation",
]

from enum import Enum
from collections.abc import Sequence

import numpy as np

from .elements import Element
from .utils import Refpts, All, refpts_iterator, _refcount

_x_axis = np.array([1.0, 0.0, 0.0])
_y_axis = np.array([0.0, 1.0, 0.0])
_z_axis = np.array([0.0, 0.0, 1.0])


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


# noinspection PyPep8Naming
def _offsets_from_translation_vector(T, ld, r3d, X_axis, Y_axis, Z_axis):
    """
    Retrieve the 3D offsets from the T1 translation vector.
    T: Translation vector (T1)
    ld: Longitudinal displacement [m]
    r3d: 3D rotation matrix
    """
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
    tD0 = T - T0
    offsets = -tD0[0] * X_axis - tD0[2] * Y_axis + (ld + tD0[5]) * Z_axis

    return offsets


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


def _ld_and_r3d_from_r_matrix(r_matrix):
    """
    Extracts the longitudinal displacement (ld) and 3D rotation matrix (r3d)
    from the R1 rotation matrix operator.
    r_matrix: 6x6 rotation matrix operator (R1).
    """
    r3d = np.eye(3)
    r3d[0, 0] = r_matrix[1, 1]
    r3d[1, 0] = r_matrix[1, 3]
    r3d[2, 0] = r_matrix[1, 4]
    r3d[0, 1] = r_matrix[3, 1]
    r3d[1, 1] = r_matrix[3, 3]
    r3d[2, 1] = r_matrix[3, 4]
    r3d[0, 2] = -r_matrix[5, 0] * r_matrix[1, 1] / r_matrix[2, 2]
    r3d[1, 2] = -r_matrix[5, 2] * r_matrix[1, 1] / r_matrix[2, 2]
    r3d[2, 2] = r_matrix[1, 1] / r_matrix[2, 2]
    ld = r_matrix[0, 1] * r_matrix[3, 3] / r_matrix[0, 0]
    return ld, r3d


class ReferencePoint(Enum):
    """Definition of the reference point for the geometric transformations."""

    CENTRE = "CENTRE"  #: Origin at the centre of the element.
    ENTRANCE = "ENTRANCE"  #: Origin at the entrance of the element.


# noinspection PyPep8Naming
def get_offsets_rotations(
    elem: Element,
    reference: ReferencePoint = ReferencePoint.CENTRE,
    *,
    RB_half: np.ndarray = None,
):
    r"""Return the offsets and rotations of a given element.

    This function returns the offsets [dx, dy, dz] and angular rotations (tilt, yaw,
    pitch) of the element.

    Args:
        elem: The beamline element.
        reference:  Transformation reference, either
                    :py:obj:`ReferencePoint.CENTRE` or
                    :py:obj:`ReferencePoint.ENTRANCE`. This must be identical to
                    the value used when the element was transformed.

    Returns:
        offsets (np.ndarray): [dx, dy, dz] array.
        tilt (float): Tilt angle (rotation about the Z-axis) in radians.
        yaw (float): Yaw angle (rotation about the Y-axis) in radians.
        pitch (float): Pitch angle (rotation about the X-axis) in radians.

    Raises:
        ValueError: If the `reference` argument is neither `ReferencePoint.CENTRE` nor
            `ReferencePoint.ENTRANCE`.

    .. Attention::

        The *reference* argument must be identical to the one used when displacing the
        element. Otherwise, the result is unpredictable.
    """
    if RB_half is None:
        elem_bending_angle = getattr(elem, "BendingAngle", 0.0)
        RB_half = _rotation([0.0, -elem_bending_angle / 2.0, 0.0])  # Eq. (12)

    try:
        T1 = elem.T1
    except AttributeError:
        T1 = np.zeros(6)
    try:
        R1 = elem.R1
    except AttributeError:
        R1 = np.eye(6)

    ld, r3d_tmp = _ld_and_r3d_from_r_matrix(R1)
    if reference is ReferencePoint.CENTRE:
        r3d = RB_half.T @ r3d_tmp @ RB_half

        X_axis = RB_half.T @ _x_axis
        Y_axis = RB_half.T @ _y_axis
        Z_axis = RB_half.T @ _z_axis

    elif reference is ReferencePoint.ENTRANCE:
        r3d = r3d_tmp

        X_axis = r3d @ _x_axis
        Y_axis = r3d @ _y_axis
        Z_axis = r3d @ _z_axis

    else:
        raise ValueError(
            "Unsupported reference, please choose either "
            "ReferencePoint.CENTRE or "
            "ReferencePoint.ENTRANCE."
        )

    offsets = _offsets_from_translation_vector(T1, ld, r3d, X_axis, Y_axis, Z_axis)
    tilt = np.arctan2(-r3d[0, 1], r3d[0, 0])
    yaw = np.arcsin(r3d[0, 2])
    pitch = np.arctan2(-r3d[1, 2], r3d[2, 2])
    return offsets, tilt, yaw, pitch


# noinspection PyPep8Naming
def transform_elem(
    elem: Element,
    reference: ReferencePoint = ReferencePoint.CENTRE,
    dx: float | None = None,
    dy: float | None = None,
    dz: float | None = None,
    tilt: float | None = None,
    pitch: float | None = None,
    yaw: float | None = None,
    *,
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

    If *relative* is :py:obj:`True`, the previous translations and angles are rebuilt
    from the *r3d* matrix and incremented by the input arguments.

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
        reference:      Transformation reference, either
                        :py:obj:`ReferencePoint.CENTRE` or
                        :py:obj:`ReferencePoint.ENTRANCE`.
        dx:             Horizontal shift [m]. Default: no change.
        dy:             Vertical shift [m]. Default: no change.
        dz:             Longitudinal shift [m]. Default: no change.
        tilt:           Tilt angle [rad]. Default: no change.
        pitch:          Pitch angle [rad]. Default: no change.
        yaw:            Yaw angle [rad]. Default: no change
        relative:       If :py:obj:`True`, the input values are added to the
          previous ones.

    .. Attention::

        When combining several transformations by using multiple calls to
        :py:func:`transform_elem`, the *reference* argument must be identical for all.
        Otherwise, the result is unpredictable.

    See Also:
        :py:func:`get_offsets_rotations`
    """
    if relative:

        def _set(ini, val):
            return ini if val is None else ini + val

    else:

        def _set(ini, val):
            return ini if val is None else val

    elem_length = getattr(elem, "Length", 0)
    elem_bending_angle = getattr(elem, "BendingAngle", 0)

    # Rotated transfer matrix linked to a bending angle
    RB = _rotation([0, -elem_bending_angle, 0])  # Eq. (12)
    RB_half = _rotation([0, -elem_bending_angle / 2, 0])  # Eq. (12)

    # Get the current transformation
    offsets0, tilt0, yaw0, pitch0 = get_offsets_rotations(
        elem, reference, RB_half=RB_half
    )

    # Apply new offsets and rotations (XYZ intrinsic order)
    offsets = np.array([_set(v0, v) for v0, v in zip(offsets0, [dx, dy, dz])])
    rotations = [_set(pitch0, pitch), _set(yaw0, yaw), _set(tilt0, tilt)]

    if reference is ReferencePoint.CENTRE:
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

    elif reference is ReferencePoint.ENTRANCE:
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


def tilt_elem(elem: Element, rots: float, relative: bool = False) -> None:
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
        relative:       If :py:obj:`True`, the rotation is added to the
          existing one

    See Also:
        :py:func:`shift_elem`
        :py:func:`.transform_elem`
    """
    transform_elem(elem, tilt=rots, relative=relative)


def shift_elem(
    elem: Element,
    dx: float | None = 0.0,
    dy: float | None = 0.0,
    dz: float | None = 0.0,
    *,
    relative: bool = False,
) -> None:
    r"""Sets the translations of an :py:class:`.Element`

    The translation vectors are stored in the :pycode:`T1` and :pycode:`T2`
    attributes.

    Parameters:
        elem:           Element to be shifted
        dx:             Horizontal translation [m]. Use :py:obj:`None` to keep
          the current value.
        dy:             Vertical translation [m]. Use :py:obj:`None` to keep
          the current value.
        dz:             Longitudinal translation [m]. Use :py:obj:`None` to keep
          the current value.
        relative:       If :py:obj:`True`, the translation is added to the
          existing one

    See Also:
        :py:func:`tilt_elem`
        :py:func:`.transform_elem`
    """
    transform_elem(elem, dx=dx, dy=dy, dz=dz, relative=relative)


def set_rotation(
    ring: Sequence[Element],
    tilts=0.0,
    pitches=0.0,
    yaws=0.0,
    *,
    refpts: Refpts = All,
    relative=False,
) -> None:
    r"""Sets the rotations of a list of elements.

    Parameters:
        ring:       Lattice description.
        tilts:      Scalar or Sequence of tilt values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        pitches:    Scalar or Sequence of pitch values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        yaws:       Scalar or Sequence of yaw values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        refpts:     Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        relative:   If :py:obj:`True`, the rotations are added to the existing ones.

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
    for el, tilt, pitch, yaw in zip(refpts_iterator(ring, refpts), tilts, pchs, yaws):
        transform_elem(el, tilt=tilt, pitch=pitch, yaw=yaw, relative=relative)


def set_tilt(
    ring: Sequence[Element], tilts, *, refpts: Refpts = All, relative=False
) -> None:
    r"""Sets the tilts of a list of elements.

    Parameters:
        ring:       Lattice description.
        tilts:      Scalar or Sequence of tilt values applied to the
          selected elements. Use :py:obj:`None` to keep the current values.
        refpts:     Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        relative:   If :py:obj:`True`, the rotation is added to the existing one.

    .. versionadded:: 0.7.0
       The *refpts* argument

    See Also:
        :py:func:`set_rotation`
        :py:func:`set_shift`
    """
    nb = _refcount(ring, refpts, endpoint=False)
    tilts = np.broadcast_to(tilts, (nb,))
    for el, tilt in zip(refpts_iterator(ring, refpts), tilts):
        transform_elem(el, tilt=tilt, relative=relative)


def set_shift(
    ring: Sequence[Element], dxs, dys, dzs=None, *, refpts: Refpts = All, relative=False
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
        transform_elem(el, dx=dx, dy=dy, dz=dz, relative=relative)


def _get_dx(elem: Element) -> float:
    """Horizontal element shift"""
    offsets, _, _, _ = get_offsets_rotations(elem, ReferencePoint.CENTRE)
    return offsets[0]


def _set_dx(elem: Element, value: float) -> None:
    transform_elem(elem, ReferencePoint.CENTRE, dx=value)


def _get_dy(elem: Element) -> float:
    """Vertical element shift"""
    offsets, _, _, _ = get_offsets_rotations(elem, ReferencePoint.CENTRE)
    return offsets[1]


def _set_dy(elem: Element, value: float) -> None:
    transform_elem(elem, ReferencePoint.CENTRE, dy=value)


def _get_dz(elem: Element) -> float:
    """Longitudinal element shift"""
    offsets, _, _, _ = get_offsets_rotations(elem, ReferencePoint.CENTRE)
    return offsets[2]


def _set_dz(elem: Element, value: float) -> None:
    transform_elem(elem, ReferencePoint.CENTRE, dz=value)


def _get_tilt(elem: Element) -> float:
    """Element tilt"""
    _, tilt, _, _ = get_offsets_rotations(elem, ReferencePoint.CENTRE)
    return tilt


def _set_tilt(elem: Element, value: float) -> None:
    transform_elem(elem, ReferencePoint.CENTRE, tilt=value)


Element.transform = transform_elem
Element.dx = property(_get_dx, _set_dx)
Element.dy = property(_get_dy, _set_dy)
Element.dz = property(_get_dz, _set_dz)
Element.tilt = property(_get_tilt, _set_tilt)
