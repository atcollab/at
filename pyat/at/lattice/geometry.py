from __future__ import annotations

__all__ = ["get_geometry"]

from collections.abc import Sequence
from math import sin, cos, atan2

import numpy as np

from .elements import Element
from .utils import Refpts, All, get_bool_index

_GEOMETRY_DTYPE = [
    ("x", np.float64),
    ("y", np.float64),
    ("z", np.float64),
    ("angle", np.float64),
]


def get_geometry(
    ring: Sequence[Element],
    refpts: Refpts = All,
    start_coordinates: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0),
):
    # noinspection PyShadowingNames
    r"""Compute the 2D ring geometry in cartesian coordinates

    Parameters:
        ring:               Lattice description.
        refpts:     Element selection key.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        start_coordinates:  *x*, *y*, *z*, *angle* at starting point. *angle* is the
          initial angle in the xy plane.

    Returns:
        geomdata:           recarray containing, x, y, z, angle.

    Example:

       >>> geomdata = get_geometry(ring)
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
            conv @= rots(elem.R1)
        if hasattr(elem, "BendingAngle"):
            ang = 0.5 * elem.BendingAngle
            rm = hkick(-ang)
            conv @= rm
            xyz += conv @ np.array([length / ang * sin(ang), 0.0, 0.0])
            conv @= rm
        else:
            xyz += conv @ np.array([length, 0.0, 0.0])
        if hasattr(elem, "R2"):
            conv @= rots(elem.R2)
        angle = atan2(conv[1, 0], conv[0, 0])
        return tuple(xyz) + (angle,)

    boolrefs = get_bool_index(ring, refpts, endpoint=True)
    x0, y0, z0, t0 = start_coordinates
    xyzc = np.array([x0, y0, z0])
    convi = hkick(t0)
    ret = [start_coordinates] + [increment(xyzc, convi, el) for el in ring]
    geomdata = np.rec.array(ret, dtype=_GEOMETRY_DTYPE)
    return geomdata[boolrefs]
