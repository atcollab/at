"""Helper functions for axis and plane descriptions"""

from __future__ import annotations

# Necessary for type aliases in python <= 3.8 :
from typing import Tuple
from typing import Union

AxisCode = Union[str, int, slice, None, type(Ellipsis)]
AxisDef = Union[AxisCode, Tuple[AxisCode, AxisCode]]

_axis_def = {
    "x": {"index": 0, "label": "x", "unit": " [m]"},
    "px": {"index": 1, "label": r"$p_x$", "unit": " [rad]"},
    "y": {"index": 2, "label": "y", "unit": " [m]"},
    "py": {"index": 3, "label": r"$p_y$", "unit": " [rad]"},
    "dp": {"index": 4, "label": r"$\delta$", "unit": ""},
    "ct": {"index": 5, "label": r"$\beta c \tau$", "unit": " [m]"},
}
for xk, xv in list(_axis_def.items()):
    xv["code"] = xk
    _axis_def[xv["index"]] = xv
    _axis_def[xk.upper()] = xv
_axis_def["delta"] = _axis_def["dp"]
_axis_def["xp"] = _axis_def["px"]  # For backward compatibility
_axis_def["yp"] = _axis_def["py"]  # For backward compatibility
_axis_def["s"] = _axis_def["ct"]
_axis_def["S"] = _axis_def["ct"]
_axis_def[None] = {"index": None, "label": "", "unit": "", "code": ":"}
_axis_def[Ellipsis] = {"index": Ellipsis, "label": "", "unit": "", "code": "..."}

_plane_def = {
    "x": {"index": 0, "label": "x", "unit": " [m]"},
    "y": {"index": 1, "label": "y", "unit": " [m]"},
    "z": {"index": 2, "label": "z", "unit": ""},
}
for xk, xv in list(_plane_def.items()):
    xv["code"] = xk
    _plane_def[xv["index"]] = xv
    _plane_def[xk.upper()] = xv
_plane_def["h"] = _plane_def["x"]
_plane_def["v"] = _plane_def["y"]
_plane_def["H"] = _plane_def["x"]
_plane_def["V"] = _plane_def["y"]
_plane_def[None] = {"index": None, "label": "", "unit": "", "code": ":"}
_plane_def[Ellipsis] = {"index": Ellipsis, "label": "", "unit": "", "code": "..."}


def _descr(dd: dict, *args: AxisDef, key: str | None = None):
    for arg in args:
        if isinstance(arg, tuple):
            for a in arg:
                yield from _descr(dd, a, key=key)
        else:
            if isinstance(arg, slice):
                descr = {"index": arg, "code": arg, "label": "", "unit": ""}
            else:
                descr = dd[arg]
            if key is None:
                yield descr
            else:
                yield descr[key]


def axis_(*axis: AxisDef, key: str | None = None):
    r"""Return axis descriptions

    Parameters:
        axis:           axis code or tuple of axis codes selecting axes in the
          standard AT coordinate system. codes for the 6 axes are:

          0, 'x', 'X'

          1, 'px', PX'

          2, 'y', 'Y'

          3, 'py', 'PY'

          4, 'dp', 'DP', 'delta'

          5, 'ct', 'CT', 's', 'S'

          :py:obj:`None` and :py:obj:`Ellipsis` select all axes

        key:            key in the coordinate descriptor dictionary,
          selecting the desired information. One of :

          'index'
            index in the standard AT coordinate vector
          'code'
            string representation
          'label'
            label for plot annotation
          'unit'
            coordinate unit
          :py:obj:`None`
            entire description dictionary

    Returns:
        descr : value or tuple[values]

    Examples:

        >>> axis_("x", "dp", key="index")
        (0, 4)

        returns the indices in the standard coordinate vector

        >>> dplabel = axis_("dp", key="label")
        >>> print(dplabel)
        $\delta$

        returns the coordinate label for plot annotation

        >>> axis_(0, "dp")
        ({'plane': 0, 'label': 'x', 'unit': ' [m]', 'code': 'x'},
         {'plane': 4, 'label': '$\\delta$', 'unit': '', 'code': 'dp'})

        returns the entire description directories

    """
    ret = tuple(_descr(_axis_def, *axis, key=key))
    if len(ret) > 1:
        return ret
    else:
        return ret[0]


def plane_(*plane: AxisDef, key: str | None = None):
    r"""Return plane descriptions

    Parameters:
        plane:          plane code or tuple of plane codes selecting planes.
          codes for the 3 planes are:

          0, 'x', 'X', 'h', 'H' for horizontal plane,

          1, 'y', 'Y', 'v', 'V' for vertical plane,

          2, 'z', 'Z'   for the longitudinal plane

          :py:obj:`None`, slice(None) and :py:obj:`Ellipsis` selects all planes
        key:            key in the plane descriptor dictionary,
          selecting the desired information. One of :

          'plane'
            plane in the optics data
          'code'
            string representation
          'label'
            label for plot annotation
          'unit'
            coordinate unit
          :py:obj:`None`
            entire description dictionary

    Returns:
        descr : value or tuple[values]

    Examples:

        >>> plane_("v", key="index")
        1

        returns the indices in the standard coordinate vector

        >>> plane_("x", "y")
        ({'plane': 0, 'label': 'x', 'unit': ' [m]', 'code': 'h'},
         {'plane': 1, 'label': 'y', 'unit': ' [m]', 'code': 'v'})

        returns the entire description directories

    """
    ret = tuple(_descr(_plane_def, *plane, key=key))
    if len(ret) > 1:
        return ret
    else:
        return ret[0]
