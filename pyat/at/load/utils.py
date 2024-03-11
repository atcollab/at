"""
Conversion utilities for creating pyat elements
"""

from __future__ import annotations

import collections
import os
import re
import sysconfig
from typing import Optional
from warnings import warn

import numpy as np

# imports necessary in' globals()' for 'eval'
from numpy import array, uint8, NaN  # noqa: F401

from at import integrators
from at.lattice import AtWarning
from at.lattice import CLASS_MAP, elements as elt
from at.lattice import Particle, Element
from at.lattice import idtable_element

_ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
_relativistic_particle = Particle()


def _particle(value) -> Particle:
    if isinstance(value, Particle):
        # Create from python: save_mat
        return value
    else:
        # Create from Matlab: load_mat
        name = value.pop("name")
        return Particle(name, **value)


def _warn(index: int, message: str, elem_dict: dict) -> None:
    name = elem_dict.get("FamName", "")
    location = f'"{name}":\n' if index is None else f'{index} ("{name}"):\n'
    warning = "".join(("In element ", location, message, f"\n{elem_dict}\n"))
    warn(AtWarning(warning), stacklevel=2)


class RingParam(elt.Element):
    """Private class for Matlab RingParam element

    :meta private:
    """

    # noinspection PyProtectedMember
    _BUILD_ATTRIBUTES = elt.Element._BUILD_ATTRIBUTES + [
        "Energy",
        "Periodicity",
    ]
    _conversions = dict(
        elt.Element._conversions, Energy=float, Periodicity=int, Particle=_particle
    )

    def __init__(
        self,
        name: str,
        energy: float,
        periodicity: int = 1,
        particle: Particle = _relativistic_particle,
        **kwargs,
    ):
        if not np.isnan(float(energy)):
            kwargs.setdefault("Energy", energy)
        kwargs.setdefault("Periodicity", periodicity)
        kwargs.setdefault("Particle", particle)
        kwargs.setdefault("PassMethod", "IdentityPass")
        super(RingParam, self).__init__(name, **kwargs)


_alias_map = {
    "rbend": elt.Dipole,
    "sbend": elt.Dipole,
    "quad": elt.Quadrupole,
    "sext": elt.Sextupole,
    "rf": elt.RFCavity,
    "bpm": elt.Monitor,
    "ap": elt.Aperture,
    "ringparam": RingParam,
    "wig": elt.Wiggler,
    "insertiondevicekickmap": idtable_element.InsertionDeviceKickMap,
    "matrix66": elt.M66,
}


# Matlab to Python class translation
_CLASS_MAP = dict((k.lower(), v) for k, v in CLASS_MAP.items())
_CLASS_MAP.update(_alias_map)

_PASS_MAP = {
    "BendLinearPass": elt.Dipole,
    "BndMPoleSymplectic4RadPass": elt.Dipole,
    "BndMPoleSymplectic4Pass": elt.Dipole,
    "QuadLinearPass": elt.Quadrupole,
    "StrMPoleSymplectic4Pass": elt.Multipole,
    "StrMPoleSymplectic4RadPass": elt.Multipole,
    "CorrectorPass": elt.Corrector,
    "CavityPass": elt.RFCavity,
    "RFCavityPass": elt.RFCavity,
    "ThinMPolePass": elt.ThinMultipole,
    "Matrix66Pass": elt.M66,
    "AperturePass": elt.Aperture,
    "IdTablePass": idtable_element.InsertionDeviceKickMap,
    "GWigSymplecticPass": elt.Wiggler,
}

# Matlab to Python attribute translation
_param_to_lattice = {
    "Energy": "energy",
    "Periodicity": "periodicity",
    "FamName": "name",
}

# Python to Matlab class translation
_matclass_map = {
    "Dipole": "Bend",
    "InsertionDeviceKickMap": "InsertionDeviceKickMap",
    "M66": "Matrix66",
}

# Python to Matlab type translation
_mattype_map = {
    int: float,
    np.ndarray: lambda attr: np.asanyarray(attr),
    Particle: lambda attr: attr.to_dict(),
}

_class_to_matfunc = {
    elt.Dipole: "atsbend",
    elt.Bend: "atsbend",
    elt.M66: "atM66",
    idtable_element.InsertionDeviceKickMap: "atinsertiondevicekickmap",
}


def hasattrs(kwargs: dict, *attributes) -> bool:
    """Checks the presence of keys in a :py:class:`dict`

    Returns :py:obj:`True` if any of the ``attributes`` is in ``kwargs``

    Args:
        kwargs:     The dictionary of keyword arguments passed to the
          Element constructor.
        attributes: A list of strings, the attribute names to be checked.

    Returns:
        found (bool):   :py:obj:`True` if the element has any of the specified
          attributes.
    """
    for attribute in attributes:
        if attribute in kwargs:
            return True
    return False


def find_class(
    elem_dict: dict, quiet: bool = False, index: Optional[int] = None
) -> type(Element):
    """Identify the class of an element from its attributes

    Args:
        elem_dict:      The dictionary of keyword arguments passed to the
                        Element constructor.
        quiet:          Suppress the warning for non-standard classes
        index:          Element index in the lattice

    Returns:
        element_class:  The guessed Class name
    """

    def low_order(key):
        polynom = np.array(elem_dict[key], dtype=np.float64).reshape(-1)
        try:
            low = np.where(polynom != 0.0)[0][0]
        except IndexError:
            low = -1
        return low

    class_name = elem_dict.get("Class", "")
    try:
        return _CLASS_MAP[class_name.lower()]
    except KeyError:
        if not quiet and class_name:
            _warn(index, f"Class '{class_name}' does not exist.", elem_dict)
        fam_name = elem_dict.get("FamName", "")
        try:
            return _CLASS_MAP[fam_name.lower()]
        except KeyError:
            pass_method = elem_dict.get("PassMethod", "")
            if not quiet and not pass_method:
                _warn(index, "No PassMethod provided.", elem_dict)
            elif not quiet and not pass_method.endswith("Pass"):
                message = (
                    f"Invalid PassMethod '{pass_method}', "
                    "provided pass methods should end in 'Pass'."
                )
                _warn(index, message, elem_dict)
            class_from_pass = _PASS_MAP.get(pass_method)
            if class_from_pass is not None:
                return class_from_pass
            else:
                length = float(elem_dict.get("Length", 0.0))
                if hasattrs(
                    elem_dict,
                    "FullGap",
                    "FringeInt1",
                    "FringeInt2",
                    "gK",
                    "EntranceAngle",
                    "ExitAngle",
                ):
                    return elt.Dipole
                elif hasattrs(
                    elem_dict,
                    "Voltage",
                    "Frequency",
                    "HarmNumber",
                    "PhaseLag",
                    "TimeLag",
                ):
                    return elt.RFCavity
                elif hasattrs(elem_dict, "Periodicity"):
                    # noinspection PyProtectedMember
                    return RingParam
                elif hasattrs(elem_dict, "Limits"):
                    return elt.Aperture
                elif hasattrs(elem_dict, "M66"):
                    return elt.M66
                elif hasattrs(elem_dict, "K"):
                    return elt.Quadrupole
                elif hasattrs(elem_dict, "PolynomB", "PolynomA"):
                    loworder = low_order("PolynomB")
                    if loworder == 1:
                        return elt.Quadrupole
                    elif loworder == 2:
                        return elt.Sextupole
                    elif loworder == 3:
                        return elt.Octupole
                    elif pass_method.startswith("StrMPoleSymplectic4") or (length > 0):
                        return elt.Multipole
                    else:
                        return elt.ThinMultipole
                elif hasattrs(elem_dict, "KickAngle"):
                    return elt.Corrector
                elif length > 0.0:
                    return elt.Drift
                elif hasattrs(elem_dict, "GCR"):
                    return elt.Monitor
                elif pass_method == "IdentityPass":
                    return elt.Marker
                else:
                    return elt.Element


def element_from_dict(
    elem_dict: dict,
    index: Optional[int] = None,
    check: bool = True,
    quiet: bool = False,
) -> Element:
    """Builds an :py:class:`.Element` from a dictionary of attributes

    Parameters:
        elem_dict:      Dictionary of element attributes
        index:          Element index
        check:          Check the compatibility of class and PassMethod
        quiet:          Suppress the warning for non-standard classes

    Returns:
        elem (Element): new :py:class:`.Element`
    """

    # noinspection PyShadowingNames
    def sanitise_class(index, cls, elem_dict):
        """Checks that the Class and PassMethod of the element are a valid
            combination. Some Classes and PassMethods are incompatible and
            would raise errors during calculation, so we send a
            warning here.

        Args:
            index:          element index
            cls:            Proposed class
            elem_dict:      The dictionary of keyword arguments passed to the
                            Element constructor.

        Raises:
            AttributeError: if the PassMethod and Class are incompatible.
        """
        class_name = cls.__name__
        pass_method = elem_dict.get("PassMethod")
        if pass_method is not None:
            pass_to_class = _PASS_MAP.get(pass_method)
            length = float(elem_dict.get("Length", 0.0))
            file_name = pass_method + _ext_suffix
            file_path = os.path.join(integrators.__path__[0], file_name)
            if not os.path.isfile(os.path.realpath(file_path)):
                message = f"PassMethod {pass_method} is missing {file_name}."
                _warn(index, message, elem_dict)
            elif (pass_method == "IdentityPass") and (length != 0.0):
                message = (
                    f"PassMethod {pass_method} is not compatible with length {length}."
                )
                _warn(index, message, elem_dict)
            elif pass_to_class is not None:
                if not issubclass(cls, pass_to_class):
                    message = (
                        f"PassMethod {pass_method} is not compatible "
                        f"with Class {class_name}."
                    )
                    _warn(index, message, elem_dict)

    cls = find_class(elem_dict, quiet=quiet, index=index)
    if check:
        sanitise_class(index, cls, elem_dict)
    # Remove mandatory attributes from the keyword arguments.
    # Create list rather than generator to ensure that elements are removed
    # from elem_dict.
    elem_args = [elem_dict.pop(attr, None) for attr in cls._BUILD_ATTRIBUTES]
    element = cls(*(arg for arg in elem_args if arg is not None), **elem_dict)
    return element


def element_from_string(elem_string: str) -> Element:
    """Builds an :py:class:`.Element` from its python :py:func:`repr` string

    Parameters:
        elem_string:    String representation of an :py:class:`.Element`

    Returns:
        elem (Element): new :py:class:`.Element`
    """
    return eval(elem_string, globals(), CLASS_MAP)


def element_from_m(line: str) -> Element:
    """Builds an :py:class:`.Element` from a line in an m-file

    Parameters:
        line:           Matlab string representation of an :py:class:`.Element`

    Returns:
        elem (Element): new :py:class:`.Element`
    """

    def argsplit(value):
        return [a.strip() for a in split_ignoring_parentheses(value, ",")]

    def makedir(mat_struct):
        """Build directory from Matlab struct arguments"""

        def pairs(it):
            while True:
                try:
                    a = next(it)
                except StopIteration:
                    break
                yield eval(a), convert(next(it))

        return dict(pairs(iter(mat_struct)))

    def makearray(mat_arr):
        """Build numpy array for Matlab array syntax"""

        def arraystr(arr):
            lns = arr.split(";")
            rr = [arraystr(v) for v in lns] if len(lns) > 1 else lns[0].split()
            return "[{0}]".format(", ".join(rr))

        return eval("array({0})".format(arraystr(mat_arr)))

    def convert(value):
        """convert Matlab syntax to numpy syntax"""
        if value.startswith("["):
            result = makearray(value[1:-1])
        elif value.startswith("struct"):
            result = makedir(argsplit(value[7:-1]))
        else:
            result = eval(value)
        return result

    left = line.index("(")
    right = line.rindex(")")
    matcls = line[:left].strip()[2:]
    cls = _CLASS_MAP[matcls]
    arguments = argsplit(line[left + 1 : right])
    ll = len(cls._BUILD_ATTRIBUTES)
    if ll < len(arguments) and arguments[ll].endswith("Pass'"):
        arguments.insert(ll, "'PassMethod'")
    args = [convert(v) for v in arguments[:ll]]
    kwargs = makedir(arguments[ll:])
    if matcls == "rbend":
        # the Matlab 'rbend' has no equivalent in PyAT. This adds parameters
        # necessary for using the python sector bend
        halfangle = 0.5 * args[2]
        kwargs.setdefault("EntranceAngle", halfangle)
        kwargs.setdefault("ExitAngle", halfangle)
    return cls(*args, **kwargs)


def element_to_dict(elem: Element) -> dict:
    """Builds the Matlab structure of an :py:class:`.Element`

    Parameters:
        elem:           :py:class:`.Element`

    Returns:
        dct (dict):     Dictionary of :py:class:`.Element` attributes
    """
    dct = dict(
        (k, _mattype_map.get(type(v), lambda attr: attr)(v)) for k, v in elem.items()
    )
    class_name = elem.__class__.__name__
    dct["Class"] = _matclass_map.get(class_name, class_name)
    return dct


def element_to_m(elem: Element) -> str:
    """Builds the Matlab-evaluable string for an :py:class:`.Element`

    Parameters:
        elem:           :py:class:`.Element`

    Returns:
        mstr (str):     Matlab string representation of the
          :py:class:`.Element` attributes
    """

    def convert(arg):
        def convert_dict(pdir):
            def scan(d):
                for k, v in d.items():
                    yield convert(k)
                    yield convert(v)

            return "struct({0})".format(", ".join(scan(pdir)))

        def convert_array(arr):
            if arr.ndim > 1:
                lns = (str(list(ln)).replace(",", "")[1:-1] for ln in arr)
                return "".join(("[", "; ".join(lns), "]"))
            elif arr.ndim > 0:
                return str(list(arr)).replace(",", "")
            else:
                return str(arr)

        if isinstance(arg, np.ndarray):
            return convert_array(arg)
        elif isinstance(arg, dict):
            return convert_dict(arg)
        elif isinstance(arg, Particle):
            return convert_dict(arg.to_dict())
        else:
            return repr(arg)

    def m_name(elclass):
        stdname = "".join(("at", elclass.__name__.lower()))
        return _class_to_matfunc.get(elclass, stdname)

    attrs = dict(elem.items())
    # noinspection PyProtectedMember
    args = [attrs.pop(k, getattr(elem, k)) for k in elem._BUILD_ATTRIBUTES]
    defelem = elem.__class__(*args)
    kwds = dict(
        (k, v)
        for k, v in attrs.items()
        if not np.array_equal(v, getattr(defelem, k, None))
    )
    argstrs = [convert(arg) for arg in args]
    if "PassMethod" in kwds:
        argstrs.append(convert(kwds.pop("PassMethod")))
    argstrs += [", ".join((repr(k), convert(v))) for k, v in kwds.items()]
    return "{0:>15}({1});...".format(m_name(elem.__class__), ", ".join(argstrs))


# Kept for compatibility but should be deprecated:


CLASS_MAPPING = dict((key, cls.__name__) for (key, cls) in _CLASS_MAP.items())

PASS_MAPPING = dict((key, cls.__name__) for (key, cls) in _PASS_MAP.items())


def find_class_name(elem_dict, quiet=False):
    """Derive the class name of an Element from its attributes"""
    return find_class(elem_dict, quiet=quiet).__name__


def split_ignoring_parentheses(string, delimiter):
    placeholder = "placeholder"
    substituted = string[:]
    matches = collections.deque(re.finditer("\\(.*?\\)", string))
    for match in matches:
        substituted = substituted.replace(match.group(), placeholder, 1)
    parts = substituted.split(delimiter)
    replaced_parts = []
    for part in parts:
        if placeholder in part:
            next_match = matches.popleft()
            part = part.replace(placeholder, next_match.group(), 1)
        replaced_parts.append(part)
    assert not matches

    return replaced_parts
