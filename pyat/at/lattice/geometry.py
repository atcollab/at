from __future__ import annotations

__all__ = ["get_geometry", "get_geometry3"]

from collections.abc import Sequence
from math import sin, cos, atan2

import numpy as np

from .elements import Element, Dipole
from .utils import Refpts, All, refpts_count, get_bool_index

_GEOMETRY_EPSIL = 1.0e-3

_GEOMETRY_DTYPE = [
    ("x", np.float64),
    ("y", np.float64),
    ("z", np.float64),
    ("angle", np.float64),
]


def get_geometry(
    ring: list[Element],
    refpts: Refpts = All,
    start_coordinates: tuple[float, float, float] = (0, 0, 0),
    centered: bool = False,
    regex: bool = False,
):
    # noinspection PyShadowingNames
    r"""Compute the 2D ring geometry in cartesian coordinates

    Parameters:
        ring:               Lattice description.
        refpts:     Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        start_coordinates:  *x*, *y*, *angle* at starting point. *x*
          and *y* are ignored if *centered* is :py:obj:`True`.
        centered:           if :py:obj:`True` the coordinates origin is the
          centre of the ring.
        regex: Use regular expression for *refpts* string matching instead of
          Unix shell-style wildcards.

    Returns:
        geomdata:           recarray containing, x, y, angle.
        radius:             machine radius at the beginning of the lattice.

            .. attention::
               This radius is different from the radius usually defined as
               :math:`C/2\pi`

    Example:

       >>> geomdata, radius = get_geometry(ring)
    """

    geom_dtype = [("x", np.float64), ("y", np.float64), ("angle", np.float64)]
    boolrefs = get_bool_index(ring, refpts, endpoint=True, regex=regex)
    nrefs = refpts_count(boolrefs, len(ring))
    geomdata = np.recarray((nrefs,), dtype=geom_dtype)
    xx = np.zeros(len(ring) + 1)
    yy = np.zeros(len(ring) + 1)
    angle = np.zeros(len(ring) + 1)
    x0, y0, t0 = start_coordinates
    x, y = 0.0, 0.0
    t = t0

    xx[0] = x
    yy[0] = y
    angle[0] = t
    for ind, el in enumerate(ring):
        ll = el.Length
        if isinstance(el, Dipole) and el.BendingAngle != 0:
            ang = 0.5 * el.BendingAngle
            ll *= np.sin(ang) / ang
        else:
            ang = 0.0
        t -= ang
        x += ll * np.cos(t)
        y += ll * np.sin(t)
        t -= ang
        xx[ind + 1] = x
        yy[ind + 1] = y
        angle[ind + 1] = t

    dff = (t + _GEOMETRY_EPSIL) % (2.0 * np.pi) - _GEOMETRY_EPSIL
    if abs(dff) < _GEOMETRY_EPSIL:
        xcenter = np.mean(xx)
        ycenter = np.mean(yy)
    elif abs(dff - np.pi) < _GEOMETRY_EPSIL:
        xcenter = 0.5 * x
        ycenter = 0.5 * y
    else:
        num = np.cos(t) * x + np.sin(t) * y
        den = np.sin(t - t0)
        xcenter = -num * np.sin(t0) / den
        ycenter = num * np.cos(t0) / den
    radius = np.sqrt(xcenter * xcenter + ycenter * ycenter)
    if centered:
        xx -= xcenter
        yy -= ycenter
    else:
        xx += x0
        yy += y0
    geomdata["x"] = xx[boolrefs]
    geomdata["y"] = yy[boolrefs]
    geomdata["angle"] = angle[boolrefs]
    return geomdata, radius


def get_geometry3(
    ring: Sequence[Element],
    refpts: Refpts = All,
    start_coordinates: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0),
):
    # noinspection PyShadowingNames
    r"""Compute the 3D ring geometry in cartesian coordinates

    Parameters:
        ring:               Lattice description.
        refpts:             Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        start_coordinates:  *x*, *y*, *z*, *angle* at starting point. *angle* is the
          initial angle in the xy plane.

    Returns:
        geomdata:           recarray containing, x, y, z, angle.

    Example:

       >>> geomdata = get_geometry3(ring)
       >>> xcoord = geomdata.x
    """

    def rots(rotmat):
        cns = rotmat[0, 0]
        sns = rotmat[0, 2]
        return np.array([[1.0, 0.0, 0.0], [0, cns, -sns], [0.0, sns, cns]])

    def hkick(ang: float):
        cns = cos(ang)
        sns = sin(ang)
        return np.array([[cns, -sns, 0.0], [sns, cns, 0.0], [0.0, 0.0, 1.0]])

    def increment(xyz, conv, elem):
        length = elem.Length
        if hasattr(elem, "R1"):
            # conv @= rots(elem.R1)
            ctemp = conv.copy()
            np.matmul(ctemp,rots(elem.R1), out=conv)
        if hasattr(elem, "BendingAngle"):
            ang = 0.5 * elem.BendingAngle
            rm = hkick(-ang)
            # conv @= rm
            ctemp = conv @ rm
            xyz += ctemp @ np.array([length / ang * sin(ang), 0.0, 0.0])
            # conv @= rm
            np.matmul(ctemp, rm, out=conv)
        else:
            xyz += conv @ np.array([length, 0.0, 0.0])
        if hasattr(elem, "R2"):
            # conv @= rots(elem.R2)
            ctemp = conv.copy()
            np.matmul(ctemp, rots(elem.R2), out=conv)
        angle = atan2(conv[1, 0], conv[0, 0])
        return tuple(xyz) + (angle,)

    boolrefs = get_bool_index(ring, refpts, endpoint=True)
    x0, y0, z0, t0 = start_coordinates
    xyzc = np.array([x0, y0, z0])
    convi = hkick(t0)
    coords = [start_coordinates] + [increment(xyzc, convi, el) for el in ring]
    geomdata = np.rec.array(coords, dtype=_GEOMETRY_DTYPE)
    return geomdata[boolrefs]
