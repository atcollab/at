import sys
from typing import Optional, Union
if sys.version_info.minor < 9:
    from typing import Tuple
else:
    Tuple = tuple

AxisCode = Union[str, int, slice, None, type(Ellipsis)]
AxisDef = Union[AxisCode, Tuple[AxisCode, AxisCode]]

_axis_def = dict(
    x=dict(index=0, label="x", unit=" [m]"),
    px=dict(index=1, label=r"$p_x$", unit=" [rad]"),
    y=dict(index=2, label="y", unit=" [m]"),
    py=dict(index=3, label=r"$p_y$", unit=" [rad]"),
    dp=dict(index=4, label=r"$\delta$", unit=""),
    ct=dict(index=5, label=r"$\beta c \tau$", unit=" [m]"),
)
for xk, xv in [it for it in _axis_def.items()]:
    xv['code'] = xk
    _axis_def[xv['index']] = xv
    _axis_def[xk.upper()] = xv
_axis_def['xp'] = _axis_def['px']
_axis_def['yp'] = _axis_def['py']
_axis_def['delta'] = _axis_def['dp']
_axis_def[None] = dict(index=slice(None), label="", unit="", code=":")
_axis_def[Ellipsis] = dict(index=Ellipsis, label="", unit="", code="...")

_plane_def = dict(
    h=dict(index=0, label="h", unit=" [m]"),
    v=dict(index=1, label="v", unit=" [m]"),
    l=dict(index=2, label="l", unit="")
)
for xk, xv in [it for it in _plane_def.items()]:
    xv['code'] = xk
    _plane_def[xv['index']] = xv
    _plane_def[xk.upper()] = xv
_plane_def['x'] = _plane_def['h']
_plane_def['y'] = _plane_def['v']
_plane_def['z'] = _plane_def['l']
_plane_def[None] = dict(index=slice(None), label="", unit="", code=":")
_plane_def[Ellipsis] = dict(index=Ellipsis, label="", unit="", code="...")

_no_def = {
    None: dict(index=slice(None), label="", unit="", code=":"),
    Ellipsis: dict(index=Ellipsis, label="", unit="", code="...")
}


def _descr(dd: dict, arg: AxisDef, key: Optional[str] = None):
    if isinstance(arg, tuple):
        return tuple(_descr(dd, a, key=key) for a in arg)
    else:
        try:
            descr = dd[arg]
        except (TypeError, KeyError):
            descr = dict(index=arg, code=arg, label="", unit="")
        if key is None:
            return descr
        else:
            return descr[key]


def axis_(axis: AxisDef, key: Optional[str] = None):
    r"""Return axis descriptions

    Parameters:
        axis:           code is either an integer in 0:6
          or a string in ['x', 'px', 'y', 'py', 'dp', 'ct']
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
        descr : requested information for each input argument.

    Examples:

        >>> axis_(('x','dp'), key='plane')
        (0, 4)

        returns the indices in the standard coordinate vector

        >>> dplabel = axis_('dp', key='label')
        >>> print(dplabel)
        $\delta$

        returns the coordinate label for plot annotation

        >>> axis_(('x','dp'))
        ({'plane': 0, 'label': 'x', 'unit': ' [m]', 'code': 'x'},
         {'plane': 4, 'label': '$\\delta$', 'unit': '', 'code': 'dp'})

        returns the entire description directories

    """
    return _descr(_axis_def, axis, key=key)


def plane_(plane: AxisDef, key: Optional[str] = None):
    r"""Return plane descriptions

    Parameters:
        plane:          plane is either an integer in 0:3 or
          a string in {'x', 'X', 'h', 'H', 'y', 'Y', 'v', 'V', 'z', 'Z'}
        key:            key in the coordinate descriptor dictionary,
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
        descr : requested information for each input argument.

    Examples:

        >>> plane_('y', key='plane')
        1

        returns the indices in the standard coordinate vector

        >>> plane_(('x','y'))
        ({'plane': 0, 'label': 'x', 'unit': ' [m]', 'code': 'x'},
         {'plane': 1, 'label': 'y', 'unit': ' [m]', 'code': 'y'})

        returns the entire description directories

    """
    return _descr(_plane_def, plane, key=key)


def optics_(param: str, axis: AxisDef, key: Optional[str] = None):
    if callable(param):
        return _descr(_no_def, axis, key=key)
    elif param in {'M', 'closed_orbit', 'dispersion', 'A', 'R'}:
        return _descr(_axis_def, axis, key=key)
    else:
        return _descr(_plane_def, axis, key=key)


def nop_(axis: AxisDef, key: Optional[str] = None):
    return _descr(_no_def, axis, key=key)
